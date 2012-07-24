#include "EMSquare2D.hpp"
#include <fstream>
#include <sstream>
#define _USE_MATH_DEFINES

EMSquare2D::EMSquare2D(
		       double L_x, double L_y, double T,
		       unsigned int n_cells, unsigned int n_steps,
		       EMSquare2DInitializer* init, mpi::communicator & world
		       )
: world(world),
  x_line(world.split(world.rank() / static_cast<unsigned int>(sqrt(world.size())))),
  y_line(world.split(world.rank() % static_cast<unsigned int>(sqrt(world.size())))),
  n_steps(n_steps),
  dx(L_x/n_cells),
  dy(L_y/n_cells),
  dt(T/n_steps),
  block_size(n_cells / floor(sqrt(world.size()))),
  E_x(new double*[block_size]),
  E_y(new double*[block_size+1]),
  B_z(new double*[block_size]),
  rhs_holder(new double[block_size]),
  B_top_bdy(new double[block_size]),
  B_bot_bdy(new double[block_size]),
  B_left_bdy(new double[block_size]),
  B_right_bdy(new double[block_size]),
  bdy_out_buffer_topright(new double[block_size]),
  bdy_out_buffer_botleft(new double[block_size]),
  Implicit_dy(y_line,
	      n_cells,
	      EMToeplitzMatrixInitializer(
					  1 + 2.0*C*C*dt*dt/((2*dy)*(2*dy)),
					  -C*C*dt*dt/((2*dy)*(2*dy)),
					  n_cells)),
  Implicit_dx(x_line,
	      n_cells,
	      EMToeplitzMatrixInitializer(
					  1 + 2.0*C*C*dt*dt/((2*dx)*(2*dx)),
					  -C*C*dt*dt/((2*dx)*(2*dx)),
					  n_cells)) {
  // Allocate storage referred to by row pointers
  // Each processor owns a block_size x block_size
  // square of cells, with the B fields at the center
  // of the cells and the E-fields at the edges
  // Note that this means each processor shares a set of
  // E-fields with its neighbors.
  E_x[0] = new double[block_size*(block_size+1)];
  E_y[0] = new double[(block_size+1)*block_size];
  B_z[0] = new double[block_size*block_size];

  // Set row pointers
  for(unsigned int i = 1; i < block_size; i++) {
    E_x[i] = E_x[i-1] + block_size + 1;
    E_y[i] = E_y[i-1] + block_size;
    B_z[i] = B_z[i-1] + block_size;
  }
  E_y[block_size] = E_y[block_size-1] + block_size;

  // Initialize Fields - there are more E-fields than B-fields
  for(unsigned int ix = 0; ix < block_size; ix++) {
    for(unsigned int iy = 0; iy < block_size; iy++) {
      E_x[ix][iy] = init->E_x(ix, iy);
      E_y[ix][iy] = init->E_y(ix, iy);
      B_z[ix][iy] = init->B_z(ix, iy);
    }
    E_x[ix][block_size] = init->E_x(ix, block_size); // Initialize top edge
  }
  for(unsigned int iy = 0; iy < block_size; iy++) {
    E_y[block_size][iy] = init->E_y(block_size, iy); // Initialize right edge
  }
}

void EMSquare2D::simulate() {
  for(unsigned int i=0; i < n_steps; i++) {
    if (i % 9 == 0) {
      std::ostringstream filename(std::ios::out);
      filename << "dump" << i / 9 << ".txt";
      this->dumpFields(filename.str());
    }
    this->TimeStep();
  }
}

void EMSquare2D::dumpFields(std::string filename) {
  int dummy;

  // If first on a row, can't start until previous row is done
  if(y_line.rank() != 0 && x_line.rank() == 0)
    world.recv(world.rank()-1, 0, dummy);

  for(unsigned int iy = 0; iy < block_size; iy++) {
    if(x_line.rank() != 0)
      x_line.recv(x_line.rank() - 1, 0, dummy);
    else if (iy != 0)
      x_line.recv(x_line.size() - 1, 0, dummy);

    std::ios_base::openmode om = (world.rank() == 0 && iy == 0) ? std::ios::out : std::ios::app;

    std::ofstream dump(filename.c_str(), om);
    for(unsigned int ix = 0; ix < block_size; ix++) {
      dump << B_z[ix][iy];
      if (ix != block_size-1 || x_line.rank() != x_line.size() - 1)
	dump << ", ";
    }
    if(x_line.rank() == x_line.size() - 1)
      dump << "\n";
    dump.flush();
    if(x_line.rank() != x_line.size() - 1)
      x_line.send(x_line.rank()+1, 0, dummy);
    else if (iy != block_size - 1)
      x_line.send(0, 0, dummy);
    else if (y_line.rank() != y_line.size() - 1)
      world.send(world.rank()+1, 0, dummy);
  }
}

void EMSquare2D::TimeStep() {
  this->implicitUpdateM();
  this->explicitUpdateP();
  this->implicitUpdateP();
  this->explicitUpdateM();
}

void EMSquare2D::printField(std::string msg) {
  if (world.rank() == 0)
    std::cout << msg;
  if(y_line.rank() == 0) {
    int dummy;
    if(x_line.rank() != 0)
      x_line.recv(x_line.rank() - 1, 0, dummy);
    for(unsigned int ix=0; ix < block_size; ix++) {
      std::cout << E_y[ix][0] << " ";
    }
    std::cout.flush();
    if (x_line.rank() == x_line.size() - 1)
      std::cout << std::endl;
    else
      x_line.send(x_line.rank() + 1, 0, dummy);
  }

  x_line.barrier();

  if (world.rank() == 0)
    std::cout << msg;
  if(y_line.rank() == 0) {
    int dummy;
    if(x_line.rank() != 0)
      x_line.recv(x_line.rank() - 1, 0, dummy);
    for(unsigned int ix=0; ix < block_size; ix++) {
      std::cout << B_z[ix][0] << " ";
    }
    std::cout.flush();
    if (x_line.rank() == x_line.size() - 1)
      std::cout << std::endl;
    else
      x_line.send(x_line.rank() + 1, 0, dummy);
  }  
}

void EMSquare2D::implicitUpdateP() {
  for(unsigned int ix=0; ix < this->block_size; ix++) {
    //Build RHS of tridiagonal equation - all necessary E values live on processor
    for(unsigned int iy=0; iy < this->block_size; iy++) {
      rhs_holder[iy] = B_z[ix][iy] + C*dt/(2*dy) * (E_x[ix][iy+1] - E_x[ix][iy]);
    }

    Implicit_dy.solve(rhs_holder);

    // Copy solution of tridiagonal equation into B array
    for(unsigned int iy=0; iy < this->block_size; iy++) {
      B_z[ix][iy] = rhs_holder[iy];
    }
  }

  // Need B values from neighboring processors to complete E update
  this->exchange_bdy_values(y_line, BDY_Y);

  for(unsigned int ix=0; ix < this->block_size; ix++) {
    // Update E - Need B on boundaries to do correctly
    for(unsigned int iy=0; iy < this->block_size + 1; iy++) {
      double B_above = (iy != block_size) ? B_z[ix][iy] : B_top_bdy[ix];
      double B_below = (iy != 0) ? B_z[ix][iy-1] : B_bot_bdy[ix];
      E_x[ix][iy] += C*dt/(2*dy) * (B_above - B_below);
    }
  }
}

void EMSquare2D::explicitUpdateP() {
  this->exchange_bdy_values(y_line, BDY_Y);
  for(unsigned int ix=0; ix < this->block_size; ix++) {
    for(unsigned int iy=0; iy < this->block_size; iy++) {
      rhs_holder[iy] = B_z[ix][iy];
      B_z[ix][iy] += C*dt/(2*dy) * (E_x[ix][iy+1] - E_x[ix][iy]);
    }
    for(unsigned int iy=0; iy < this->block_size+1; iy++) {
      double B_above = (iy != block_size) ? rhs_holder[iy] : B_top_bdy[ix];
      double B_below = (iy != 0) ? rhs_holder[iy-1] : B_bot_bdy[ix];
      E_x[ix][iy] += C*dt/(2*dy) * (B_above - B_below);
    }
  }
}

void EMSquare2D::implicitUpdateM() {
  for(unsigned int iy=0; iy < this->block_size; iy++) {
    for(unsigned int ix=0; ix < this->block_size; ix++) {
      rhs_holder[ix] = B_z[ix][iy] - C*dt/(2*dy) * (E_y[ix+1][iy] - E_y[ix][iy]);
    }
    
    Implicit_dx.solve(rhs_holder);

    // Copy solution of tridiagonal equation into B array
    for(unsigned int ix=0; ix < this->block_size; ix++) {
      B_z[ix][iy] = rhs_holder[ix];
    }
  }

  this->exchange_bdy_values(x_line, BDY_X);
  for(unsigned int iy=0; iy < this->block_size; iy++) {
    // Update E
    for(unsigned int ix=0; ix < this->block_size+1; ix++) {
      double B_left = (ix != 0) ? B_z[ix-1][iy] : B_left_bdy[iy];
      double B_right = (ix != block_size) ? B_z[ix][iy] : B_right_bdy[iy];
      E_y[ix][iy] -= C*dt/(2*dy) * (B_right - B_left);
    }
  }
}

void EMSquare2D::explicitUpdateM() {
  this->exchange_bdy_values(x_line, BDY_X);
  for(unsigned int iy=0; iy < this->block_size; iy++) {
    for(unsigned int ix=0; ix < this->block_size; ix++) {
      rhs_holder[ix] = B_z[ix][iy];
      B_z[ix][iy] -= C*dt/(2*dx) * (E_y[ix+1][iy] - E_y[ix][iy]);
    }
    for(unsigned int ix=0; ix < this->block_size+1; ix++) {
      double B_left = (ix != 0) ? rhs_holder[ix-1] : B_left_bdy[iy];
      double B_right = (ix != block_size) ? rhs_holder[ix] : B_right_bdy[iy];
      E_y[ix][iy] -= C*dt/(2*dx) * (B_right - B_left);
    }
  }
}

//FIXME - use enum instead of bool
void pass_vector(mpi::communicator comm,
		      double* outgoing, double* incoming, unsigned int n,
		      bool upwards) {
  int p = comm.rank();
  int n_p = comm.size();
  int dest = p + (upwards ? 1 : -1);
  int src = p + (upwards ? -1 : 1);


  if(p % 2 == 0) {
    if(src >= 0 && src < n_p)
      comm.recv(src, 0, incoming, n);
    if(dest < n_p && dest >= 0)
      comm.send(dest, 0, outgoing, n);
  } else {
    if (dest < n_p && dest >= 0)
      comm.send(dest, 0, outgoing, n);
    if(src >= 0 && src < n_p)
      comm.recv(src, 0, incoming, n);
  }
}

void EMSquare2D::exchange_bdy_values(mpi::communicator comm, bdy_dir dir) {
  if(dir == BDY_X) {
    for(unsigned int iy=0; iy < this->block_size; iy++) {
      bdy_out_buffer_topright[iy] = B_z[block_size-1][iy];
      bdy_out_buffer_botleft[iy] = B_z[0][iy];
    }
    pass_vector(comm, bdy_out_buffer_topright, B_left_bdy, block_size, true);
    pass_vector(comm, bdy_out_buffer_botleft, B_right_bdy, block_size, false);
    if(comm.rank() == 0)
      std::copy(bdy_out_buffer_botleft, bdy_out_buffer_botleft + block_size, B_left_bdy);
    if(comm.rank() == comm.size() - 1)
      std::copy(bdy_out_buffer_topright, bdy_out_buffer_topright + block_size, B_right_bdy);
  } else {
    for(unsigned int ix=0; ix < this->block_size; ix++) {
      bdy_out_buffer_topright[ix] = B_z[ix][block_size-1];
      bdy_out_buffer_botleft[ix] = B_z[ix][0];
    }
    pass_vector(comm, bdy_out_buffer_topright, B_bot_bdy, block_size, true);
    pass_vector(comm, bdy_out_buffer_botleft, B_top_bdy, block_size, false);
    if(comm.rank() == 0)
      std::copy(bdy_out_buffer_botleft, bdy_out_buffer_botleft + block_size, B_bot_bdy);
    if(comm.rank() == comm.size() - 1)
      std::copy(bdy_out_buffer_topright, bdy_out_buffer_topright + block_size, B_top_bdy);
  }
}


EMSquare2DInitializer::EMSquare2DInitializer(double dx, double dy, double L_x, double L_y,
					     unsigned int block_size, const mpi::communicator & comm)
: dx(dx), dy(dy),
  L_x(L_x), L_y(L_y) {
  unsigned int procs_per_edge = floor(sqrt(comm.size()));
  this->x_offset = block_size*(comm.rank() % procs_per_edge)*dx;
  this->y_offset = block_size*(comm.rank() / procs_per_edge)*dy;
}

//FIXME - Really needs to be different for E and B
double EMSquare2DInitializer::x_for_B(unsigned int i) {
  return this->x_offset + (i + .5)*dx;
}

double EMSquare2DInitializer::y_for_B(unsigned int j) {
  return this->y_offset + (j + .5)*dy;
}

double EMSquare2DInitializer::x_for_E(unsigned int i) {
  return this->x_offset + i*dx;
}

double EMSquare2DInitializer::y_for_E(unsigned int j) {
  return this->y_offset + j*dy;
}
