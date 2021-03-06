#include "Simulation.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <boost/timer.hpp>
#include <sys/time.h>
#define _USE_MATH_DEFINES

BoundaryLocation determineBoundary(mpi::communicator& world) {
  if(world.rank() == 0 && world.rank() == world.size()-1)
    return DOUBLE_BDY;
  else if(world.rank() == 0)
    return LOWER_BDY;
  else if(world.rank() == world.size() - 1)
    return UPPER_BDY;
  return NO_BDY;
}

Simulation::Simulation(
		       double L_x, double L_y, double T,
		       unsigned int n_cells, unsigned int n_steps,
		       unsigned int procs_x, unsigned int procs_y,
		       unsigned int block_size,
		       std::string& dump_dir,
		       SimulationInitializer* init, mpi::communicator & world
		       )
: world(world),
  xLine(world.split(world.rank() / procs_x)),
  yLine(world.split(world.rank() % procs_x)),
  nSteps(n_steps),
  dx(L_x/n_cells),
  dy(L_y/n_cells),
  dt(T/n_steps),
  blockSize(block_size),
  preFactorX(LIGHTSPEED*dt/(2*dx)),
  preFactorY(LIGHTSPEED*dt/(2*dy)),
  E_x(new double*[blockSize]),
  E_y(new double*[blockSize+1]),
  B_z(new double*[blockSize]),
  rhss(new double*[blockSize]),
  BTopBdy(new double[blockSize]),
  BBotBdy(new double[blockSize]),
  BLeftBdy(new double[blockSize]),
  BRightBdy(new double[blockSize]),
  BdyOutBufferTopRight(new double[blockSize]),
  BdyOutBufferBotLeft(new double[blockSize]),
  dumpDir(dump_dir) {

  VacuumMatrixInitializer mat_init_x = VacuumMatrixInitializer(dx, dt, blockSize, determineBoundary(xLine));
  VacuumMatrixInitializer mat_init_y = VacuumMatrixInitializer(dy, dt, blockSize, determineBoundary(yLine));
  VacuumCouplingInitializer coupling_init_x = VacuumCouplingInitializer(& mat_init_x, blockSize, xLine);
  VacuumCouplingInitializer coupling_init_y = VacuumCouplingInitializer(& mat_init_y, blockSize, yLine);


  std::vector<AbstractMatrixInitializer*> mat_inits_x(blockSize, & mat_init_x);
  std::vector<AbstractMatrixInitializer*> mat_inits_y(blockSize, & mat_init_y);
  std::vector<AbstractCouplingInitializer*> coupling_inits_x(blockSize, & coupling_init_x);
  std::vector<AbstractCouplingInitializer*> coupling_inits_y(blockSize, & coupling_init_y);


  xUpdateRHSs = init->initCollection(mat_inits_x, coupling_inits_x, blockSize, xLine);
  yUpdateRHSs = init->initCollection(mat_inits_y, coupling_inits_y, blockSize, yLine);

  this->allocate_fields(init);
}

void Simulation::allocate_fields(SimulationInitializer* init) {
  // Allocate storage referred to by row pointers
  // Each processor owns a blockSize x blockSize
  // square of cells, with the B fields at the center
  // of the cells and the E-fields at the edges
  // Note that this means each processor shares a set of
  // E-fields with its neighbors.
  E_x[0] = new double[(blockSize+1)*blockSize];
  E_y[0] = new double[blockSize*(blockSize+1)];
  B_z[0] = new double[blockSize*blockSize];
  rhss[0] = new double[blockSize*blockSize];

  // Set row pointers
  for(unsigned int i = 1; i < blockSize; i++) {
    E_x[i] = E_x[i-1] + blockSize + 1;
    E_y[i] = E_y[i-1] + blockSize;
    B_z[i] = B_z[i-1] + blockSize;
    rhss[i] = rhss[i-1] + blockSize;
  }
  E_y[blockSize] = E_y[blockSize-1] + blockSize;

  // Initialize Fields - there are more E-fields than B-fields
  for(unsigned int ix = 0; ix < blockSize; ix++) {
    for(unsigned int iy = 0; iy < blockSize; iy++) {
      E_x[ix][iy] = init->E_x(ix, iy);
      E_y[ix][iy] = init->E_y(ix, iy);
      B_z[ix][iy] = init->B_z(ix, iy);
    }
    E_x[ix][blockSize] = init->E_x(ix, blockSize); // Initialize top edge
  }
  for(unsigned int iy = 0; iy < blockSize; iy++) {
    E_y[blockSize][iy] = init->E_y(blockSize, iy); // Initialize right edge
  }  
}

void Simulation::simulate(bool dump, unsigned int dump_periodicity, unsigned int total_dumps) {
  unsigned int n_dumps = 0;
  timeval t1, t2;
  gettimeofday(& t1, NULL);

  std::cout << std::setprecision(10);
  for(unsigned int i=0; i <= nSteps; i++) {
    if (i % dump_periodicity == 0 && n_dumps < total_dumps && dump) {
      std::ostringstream filename(std::ios::out);
      filename << dumpDir << "/dump" << i / dump_periodicity << "dt" << 10.0/nSteps << "dx" << dx << ".txt";
      this->dumpFields(filename.str());
      n_dumps++;
    }
    if (i % 1000 == 0 && world.rank() == 0) {
      gettimeofday(& t2, NULL);
      if(i > 0)
	std::cout << 1000000*(t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) << std::endl;
      t1=t2;
    }
    this->TimeStep();
  }
}

void Simulation::dumpFields(std::string filename) {
  int dummy;

  // If first on a row, can't start until previous row is done
  if(yLine.rank() != 0 && xLine.rank() == 0)
    world.recv(world.rank()-1, 0, dummy);

  for(unsigned int iy = 0; iy < blockSize; iy++) {
    if(xLine.rank() != 0)
      xLine.recv(xLine.rank() - 1, 0, dummy);
    else if (iy != 0)
      xLine.recv(xLine.size() - 1, 0, dummy);

    std::ios_base::openmode om = (world.rank() == 0 && iy == 0) ? std::ios::out : std::ios::app;

    std::ofstream dump(filename.c_str(), om);
    dump << std::setprecision(20);
    
    for(unsigned int ix = 0; ix < blockSize; ix++) {
      dump << B_z[ix][iy];
      if (ix != blockSize-1 || xLine.rank() != xLine.size() - 1)
	dump << ", ";
    }
    if(xLine.rank() == xLine.size() - 1)
      dump << "\n";
    dump.flush();
    if(xLine.rank() != xLine.size() - 1)
      xLine.send(xLine.rank()+1, 0, dummy);
    else if (iy != blockSize - 1)
      xLine.send(0, 0, dummy);
    else if (yLine.rank() != yLine.size() - 1)
      world.send(world.rank()+1, 0, dummy);
  }
}

void Simulation::TimeStep() {
#ifndef YEE
#ifndef NO_IMPLICIT_SOLVE
  this->implicitUpdateM();
#endif
#ifndef NO_EXPLICIT_SOLVE
  this->explicitUpdateP();
#endif
#ifndef NO_IMPLICIT_SOLVE
  this->implicitUpdateP();
#endif
#ifndef NO_EXPLICIT_SOLVE
  this->explicitUpdateM();
#endif
#else
  this->yeeUpdate();
#endif
}

void Simulation::printField(std::string msg) {
  if (world.rank() == 0)
    std::cout << msg;
  if(yLine.rank() == 0) {
    int dummy;
    if(xLine.rank() != 0)
      xLine.recv(xLine.rank() - 1, 0, dummy);
    for(unsigned int ix=0; ix < blockSize; ix++) {
      std::cout << E_y[ix][0] << " ";
    }
    std::cout.flush();
    if (xLine.rank() == xLine.size() - 1)
      std::cout << std::endl;
    else
      xLine.send(xLine.rank() + 1, 0, dummy);
  }

  xLine.barrier();

  if (world.rank() == 0)
    std::cout << msg;
  if(yLine.rank() == 0) {
    int dummy;
    if(xLine.rank() != 0)
      xLine.recv(xLine.rank() - 1, 0, dummy);
    for(unsigned int ix=0; ix < blockSize; ix++) {
      std::cout << B_z[ix][0] << " ";
    }
    std::cout.flush();
    if (xLine.rank() == xLine.size() - 1)
      std::cout << std::endl;
    else
      xLine.send(xLine.rank() + 1, 0, dummy);
  }  
}

void Simulation::implicitUpdateM() {
#ifndef TOTAL_REDUCED_ONLY
#ifndef COMMUNICATION_ONLY  
  for(unsigned int iy=0; iy < this->blockSize; iy++) {
    for(unsigned int ix=0; ix < this->blockSize; ix++) {
      rhss[iy][ix] = B_z[ix][iy] - LIGHTSPEED*dt/(2*dy) * (E_y[ix+1][iy] - E_y[ix][iy]);
    }
  }
#endif
#endif

  xUpdateRHSs->doLines(rhss);

#ifndef TOTAL_REDUCED_ONLY
#ifndef COMMUNICATION_ONLY
  for(unsigned int iy=0; iy < this->blockSize; iy++) {
    for(unsigned int ix=0; ix < this->blockSize; ix++) {
      B_z[ix][iy] = rhss[iy][ix];
    }
  }
#endif  

  this->implicitMSubstituteB();
#endif
}

void Simulation::implicitUpdateP() {
#ifndef TOTAL_REDUCED_ONLY
#ifndef COMMUNICATION_ONLY
  for(unsigned int ix=0; ix < this->blockSize; ix++) {
    for(unsigned int iy=0; iy < this->blockSize; iy++) {
      rhss[ix][iy] = B_z[ix][iy] + LIGHTSPEED*dt/(2*dy) * (E_x[ix][iy+1] - E_x[ix][iy]);
    }
  }
#endif
#endif

  yUpdateRHSs->doLines(rhss);

#ifndef TOTAL_REDUCED_ONLY
#ifndef COMMUNICATION_ONLY
  std::copy(rhss[0],  rhss[0] + blockSize*blockSize, B_z[0]);
#endif

  this->implicitPSubstituteB();
#endif
}

void Simulation::implicitPSubstituteB() {
  // Need B values from neighboring processors to complete E update
#ifndef NO_NEAREST_NEIGHBOR
  this->exchangeBdyValues(yLine, BDY_Y);
#endif

#ifndef COMMUNICATION_ONLY
  for(unsigned int ix=0; ix < this->blockSize; ix++) {
    // Update E - Need B on boundaries to do correctly
    for(unsigned int iy=0; iy < this->blockSize + 1; iy++) {
      double B_above = (iy != blockSize) ? B_z[ix][iy] : BTopBdy[ix];
      double B_below = (iy != 0) ? B_z[ix][iy-1] : BBotBdy[ix];
      E_x[ix][iy] += LIGHTSPEED*dt/(2*dy) * (B_above - B_below);
    }
  }
#endif
}

void Simulation::yeeUpdate() {  
#ifndef YEE_NO_COMM
  this->exchangeBdyValues(xLine, BDY_X);
  this->exchangeBdyValues(yLine, BDY_Y);
#endif

  for(unsigned int ix=0; ix < this->blockSize; ix++) {    
    for(unsigned int iy=0; iy < this->blockSize; iy++) {
      B_z[ix][iy] += preFactorX * ((E_x[ix][iy+1] - E_x[ix][iy]) - (E_y[ix+1][iy] - E_y[ix][iy]));
    }
  }
  for(unsigned int ix=0; ix < blockSize+1; ix++) {
    for(unsigned int iy=0; iy < this->blockSize+1; iy++) {
      if(iy < blockSize) {
	double B_left = (ix != 0) ? B_z[ix-1][iy] : BLeftBdy[iy];
	double B_right = (ix != blockSize) ? B_z[ix][iy] : BRightBdy[iy];
	E_y[ix][iy] -= preFactorX * (B_right - B_left);
      }
      if(ix < blockSize) {
	double B_above = (iy != blockSize) ? B_z[ix][iy] : BTopBdy[ix];
	double B_below = (iy != 0) ? B_z[ix][iy-1] : BBotBdy[ix];
	E_x[ix][iy] += preFactorY * (B_above - B_below);
      }
    }
  }
}

void Simulation::explicitUpdateP() {
#ifndef NO_NEAREST_NEIGHBOR
  this->exchangeBdyValues(yLine, BDY_Y);
#endif
#ifndef COMMUNICATION_ONLY
  for(unsigned int ix=0; ix < this->blockSize; ix++) {    
    for(unsigned int iy=0; iy < this->blockSize; iy++) {
      // rhss[0][iy] = B_z[ix][iy];
      rhss[ix][iy] = B_z[ix][iy];
      B_z[ix][iy] += LIGHTSPEED*dt/(2*dy) * (E_x[ix][iy+1] - E_x[ix][iy]);
    }
    for(unsigned int iy=0; iy < this->blockSize+1; iy++) {
      // double B_above = (iy != blockSize) ? rhss[0][iy] : BTopBdy[ix];
      // double B_below = (iy != 0) ? rhss[0][iy-1] : BBotBdy[ix];
      double B_above = (iy != blockSize) ? rhss[ix][iy] : BTopBdy[ix];
      double B_below = (iy != 0) ? rhss[ix][iy-1] : BBotBdy[ix];
      E_x[ix][iy] += LIGHTSPEED*dt/(2*dy) * (B_above - B_below);
    }
  }
#endif
}
void Simulation::implicitMSubstituteB() {
#ifndef NO_NEAREST_NEIGHBOR
  this->exchangeBdyValues(xLine, BDY_X);
#endif
#ifndef COMMUNICATION_ONLY
  for(unsigned int iy=0; iy < this->blockSize; iy++) {
    // Update E
    for(unsigned int ix=0; ix < this->blockSize+1; ix++) {
      double B_left = (ix != 0) ? B_z[ix-1][iy] : BLeftBdy[iy];
      double B_right = (ix != blockSize) ? B_z[ix][iy] : BRightBdy[iy];
      E_y[ix][iy] -= LIGHTSPEED*dt/(2*dy) * (B_right - B_left);
    }
  }
#endif
}

void Simulation::explicitUpdateM() {
#ifndef NO_NEAREST_NEIGHBOR
  this->exchangeBdyValues(xLine, BDY_X);
#endif
#ifndef COMMUNICATION_ONLY
  for(unsigned int iy=0; iy < this->blockSize; iy++) {
    for(unsigned int ix=0; ix < this->blockSize; ix++) {
      // rhss[0][ix] = B_z[ix][iy];
      rhss[iy][ix] = B_z[ix][iy];
      B_z[ix][iy] -= LIGHTSPEED*dt/(2*dx) * (E_y[ix+1][iy] - E_y[ix][iy]);
    }
    for(unsigned int ix=0; ix < this->blockSize+1; ix++) {
      // double B_left = (ix != 0) ? rhss[0][ix-1] : BLeftBdy[iy];
      // double B_right = (ix != blockSize) ? rhss[0][ix] : BRightBdy[iy];
      double B_left = (ix != 0) ? rhss[iy][ix-1] : BLeftBdy[iy];
      double B_right = (ix != blockSize) ? rhss[iy][ix] : BRightBdy[iy];
      E_y[ix][iy] -= LIGHTSPEED*dt/(2*dx) * (B_right - B_left);
    }
  }
#endif
}

//FIXME - use enum instead of bool
void pass_vector(mpi::communicator comm,
		      double* outgoing, double* incoming, unsigned int n,
		      bool upwards) {
  int p = comm.rank();
  int n_p = comm.size();
  int partner = p + (upwards ? 1 : -1);
  // int src = p + (upwards ? -1 : 1);


  if(p % 2 == 0) {
    if(partner >= 0 && partner < n_p) {
      comm.recv(partner, 0, incoming, n);
      comm.send(partner, 0, outgoing, n);
    }
  } else {
    if (partner >= 0 && partner < n_p) {
      comm.send(partner, 0, outgoing, n);
      comm.recv(partner, 0, incoming, n);
    }
  }
}

void Simulation::exchangeBdyValues(mpi::communicator comm, bdyDir dir) {
  unsigned int p = world.rank();
  if(dir == BDY_X) {
    for(unsigned int iy=0; iy < this->blockSize; iy++) {
      BdyOutBufferTopRight[iy] = B_z[blockSize-1][iy];
      BdyOutBufferBotLeft[iy] = B_z[0][iy];
    }
    if(p % 2 == 0) {
      // std::cout << "First pass_vector on even rank" << p << std::endl;
      pass_vector(comm, BdyOutBufferTopRight, BLeftBdy, blockSize, true);
      // std::cout << "Second pass_vector on even rank" << p << std::endl;
      pass_vector(comm, BdyOutBufferBotLeft, BRightBdy, blockSize, false);
      // std::cout << "even rank done" << p <<std::endl;
    } else {
      // std::cout << "First pass_vector on odd rank" << p << std::endl;
      pass_vector(comm, BdyOutBufferBotLeft, BRightBdy, blockSize, false);
      // std::cout << "Second pass_vector on odd rank" << p << std::endl;
      pass_vector(comm, BdyOutBufferTopRight, BLeftBdy, blockSize, true);
      // std::cout << "odd rank done" << p <<std::endl;
    }
    if(comm.rank() == 0)
      std::copy(BdyOutBufferBotLeft, BdyOutBufferBotLeft + blockSize, BLeftBdy);
    if(comm.rank() == comm.size() - 1)
      std::copy(BdyOutBufferTopRight, BdyOutBufferTopRight + blockSize, BRightBdy);
  } else {
    for(unsigned int ix=0; ix < this->blockSize; ix++) {
      BdyOutBufferTopRight[ix] = B_z[ix][blockSize-1];
      BdyOutBufferBotLeft[ix] = B_z[ix][0];
    }
    if (p % 2 == 0) {
      // std::cout << "First pass_vector on even rank" << p << std::endl;
      pass_vector(comm, BdyOutBufferTopRight, BBotBdy, blockSize, true);
      // std::cout << "Second pass_vector on even rank" << p << std::endl;
      pass_vector(comm, BdyOutBufferBotLeft, BTopBdy, blockSize, false);
      // std::cout << "even rank done" << p <<std::endl;
    } else {
      // std::cout << "First pass_vector on odd rank" << p << std::endl;
      pass_vector(comm, BdyOutBufferBotLeft, BTopBdy, blockSize, false);
      // std::cout << "Second pass_vector on odd rank" << p << std::endl;
      pass_vector(comm, BdyOutBufferTopRight, BBotBdy, blockSize, true);
      // std::cout << "odd rank done" << p <<std::endl;
    }
    if(comm.rank() == 0)
      std::copy(BdyOutBufferBotLeft, BdyOutBufferBotLeft + blockSize, BBotBdy);
    if(comm.rank() == comm.size() - 1)
      std::copy(BdyOutBufferTopRight, BdyOutBufferTopRight + blockSize, BTopBdy);
  }
}

SimulationInitializer::SimulationInitializer(double dx, double dy, double L_x, double L_y,
					     unsigned int x_procs, unsigned int y_procs,
					     unsigned int blockSize, const mpi::communicator & comm)
: dx(dx), dy(dy),
  L_x(L_x), L_y(L_y),
  x_procs(x_procs), y_procs(y_procs) {
  unsigned int procs_per_edge = floor(sqrt(comm.size()));
  this->x_offset = blockSize*(comm.rank() % x_procs)*dx;
  this->y_offset = blockSize*(comm.rank() / x_procs)*dy;
}

double SimulationInitializer::x_for_B(unsigned int i) {
  return this->x_offset + (i + .5)*dx;
}

double SimulationInitializer::y_for_B(unsigned int j) {
  return this->y_offset + (j + .5)*dy;
}

double SimulationInitializer::x_for_E(unsigned int i) {
  return this->x_offset + i*dx;
}

double SimulationInitializer::y_for_E(unsigned int j) {
  return this->y_offset + j*dy;
}
