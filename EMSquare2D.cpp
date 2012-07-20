#include "EMSquare2D.hpp"
#include <fstream>
#include<sstream>
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
  E_y(new double*[block_size]),
  B_z(new double*[block_size]),
  rhs_holder(new double[block_size]),
  boundary_E_in(new double[block_size]),
  boundary_E_out(new double[block_size]),
  boundary_B_in(new double[block_size]),
  boundary_B_out(new double[block_size]),
  Implicit_dy(x_line,
			  n_cells,
			  ToeplitzMatrixInitializer(
										  1 + 2.0*C*C*dt*dt/((2*dy)*(2*dy)),
										  -C*C*dt*dt/((2*dy)*(2*dy)))),
  Implicit_dx(y_line,
			  n_cells,
			  ToeplitzMatrixInitializer(
										  1 + 2.0*C*C*dt*dt/((2*dy)*(2*dy)),
										  -C*C*dt*dt/((2*dy)*(2*dy)))) {
  E_x[0] = new double[block_size*block_size];
  E_y[0] = new double[block_size*block_size];
  B_z[0] = new double[block_size*block_size];
  for(unsigned int i = 1; i < block_size; i++) {
	E_x[i] = E_x[i-1] + block_size;
	E_y[i] = E_y[i-1] + block_size;
	B_z[i] = B_z[i-1] + block_size;
  }

  for(unsigned int ix = 0; ix < block_size; ix++) {
	for(unsigned int iy = 0; iy < block_size; iy++) {
	  E_x[ix][iy] = init->E_x(ix, iy);
	  E_y[ix][iy] = init->E_y(ix, iy);
	  B_z[ix][iy] = init->B_z(ix, iy);
	}
  }
}

void EMSquare2D::simulate() {
  for(unsigned int i=0; i < n_steps; i++) {
	if (i % 1 == 0) {
	  std::ostringstream filename(std::ios::out);
	  filename << "herp" << i / 1 << ".txt";
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
	dump << "<<" << x_line.rank() << ">>";
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
  this->printField("Inital: ");
  this->implicitUpdateM();
  world.barrier();
  this->printField("Imp M: ");
  this->explicitUpdateP();
  world.barrier();
  this->printField("Exp P: ");
  this->implicitUpdateP();
  world.barrier();
  this->printField("Imp P: ");
  world.barrier();
  this->explicitUpdateM();
  world.barrier();
  this->printField("Exp M: ");
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
  this->exchange_bdy_values(y_line, BDY_Y);
  for(unsigned int ix=0; ix < this->block_size; ix++) {
	//Build RHS of tridiagonal equation - need one B value from neighboring processor
	for(unsigned int iy=1; iy < this->block_size; iy++) {
	  rhs_holder[iy] = E_x[ix][iy] - C*dt/(2*dy) * (B_z[ix][iy] - B_z[ix][iy-1]);
	}
	rhs_holder[0] = E_x[ix][0] - C*dt/(2*dy) * (B_z[ix][0] - this->boundary_B_in[ix]);

	Implicit_dy.solve(rhs_holder);

	// Copy solution of tridiagonal equation into E array
	for(unsigned int iy=0; iy < this->block_size; iy++) {
	  E_x[ix][iy] = rhs_holder[iy];
	}
  }

  this->exchange_bdy_values(y_line, BDY_Y);
  for(unsigned int ix=0; ix < this->block_size; ix++) {
  	// Update B
	for(unsigned int iy=0; iy < this->block_size - 1; iy++) {
	  B_z[ix][iy] += C*dt/(2*dy) * (E_x[ix][iy+1] - E_x[ix][iy]);
	}
	B_z[ix][block_size-1] += C*dt/(2*dy) *(this->boundary_E_in[ix] - E_x[ix][block_size-1]);
  }
}

void EMSquare2D::explicitUpdateP() {
  this->exchange_bdy_values(y_line, BDY_Y);
  for(unsigned int ix=0; ix < this->block_size; ix++) {
	rhs_holder[0]=E_x[ix][0];
	E_x[ix][0]+=C*dt/(2*dy) * (B_z[ix][0] - boundary_B_in[ix]);
	for(unsigned int iy=1; iy < this->block_size; iy++) {
	  rhs_holder[iy] = E_x[ix][iy];
	  E_x[ix][iy] += C*dt/(2*dy) * (B_z[ix][iy] - B_z[ix][iy-1]);
	}
	for(unsigned int iy=0; iy < this->block_size-1; iy++) {
	  B_z[ix][iy] += C*dt/(2*dy) * (rhs_holder[iy+1] - rhs_holder[iy]);
	}
	B_z[ix][this->block_size-1] += C*dt/(2*dy)*(boundary_E_in[ix] - E_x[ix][this->block_size-1]);
  }
}

void EMSquare2D::implicitUpdateM() {

  // if(world.rank() == 0)
  // 	for(unsigned int iy=0; iy < block_size; iy++)
  // 	  std::cout << boundary_B_in[iy] << " " << (iy == block_size - 1 ? "\n" : "");

  this->exchange_bdy_values(x_line, BDY_X);

  // if(world.rank() == 0)
  // 	for(unsigned int iy=0; iy < block_size; iy++)
  // 	  std::cout << boundary_B_in[iy] << " " << (iy == block_size - 1 ? "\n" : "");

  for(unsigned int iy=0; iy < this->block_size; iy++) {
	// if(x_line.rank() == 0 && iy == 5)
	//   std::cout << "BOUNDARY VALUE: " << this->boundary_B_in[iy] << std::endl;
	//Build RHS of tridiagonal equation - need one B value from neighboring processor
	for(unsigned int ix=1; ix < this->block_size; ix++) {
	  rhs_holder[ix] = E_y[ix][iy] + C*dt/(2*dy) * (B_z[ix][iy] - B_z[ix-1][iy]);
	}
	rhs_holder[0] = E_y[0][iy] + C*dt/(2*dy) * (B_z[0][iy] - this->boundary_B_in[iy]);

	// if(y_line.rank() == 0 && iy == 0) {
	//   int dummy;
	//   if(x_line.rank() != 0)
	// 	x_line.recv(x_line.rank()-1, 0, dummy);
	//   else
	// 	std::cout << "pre-solve: ";
	//   for(int ix=0; ix < this->block_size; ix++)
	// 	std::cout << rhs_holder[ix] << " ";
	//   std::cout.flush();
	//   if(x_line.rank() != x_line.size()-1)
	// 	x_line.send(x_line.rank()+1, 0, dummy);
	//   else {
	// 	std::cout << std::endl;
	// 	x_line.send(0, 0, dummy);
	//   }
	// }
	
	Implicit_dx.solve(rhs_holder);

	// if(y_line.rank() == 0 && iy == 0) {
	//   int dummy;
	//   if(x_line.rank() != 0)
	// 	x_line.recv(x_line.rank()-1, 0, dummy);
	//   else {
	// 	x_line.recv(x_line.size()-1, 0, dummy);
	// 	std::cout << "post-solve: ";
	//   }
	//   for(int ix=0; ix <  this->block_size; ix++)
	// 	std::cout << rhs_holder[ix] << " ";
	//   std::cout.flush();
	//   if(x_line.rank() != x_line.size()-1)
	// 	x_line.send(x_line.rank()+1, 0, dummy);
	//   else
	// 	std::cout << std::endl;
	// }

	// Copy solution of tridiagonal equation into E array
	for(unsigned int ix=0; ix < this->block_size; ix++) {
	  E_y[ix][iy] = rhs_holder[ix];
	}
  }

  this->exchange_bdy_values(x_line, BDY_X);
  for(unsigned int iy=0; iy < this->block_size; iy++) {
  	// Update B
	for(unsigned int ix=0; ix < this->block_size - 1; ix++) {
	  B_z[ix][iy] += C*dt/(2*dy) * (E_y[ix+1][iy] - E_y[ix][iy]);
	}
	B_z[block_size-1][iy] += C*dt/(2*dy) *(this->boundary_E_in[iy] - E_y[block_size-1][iy]);
  }
}

void EMSquare2D::explicitUpdateM() {
  this->exchange_bdy_values(x_line, BDY_X);
  for(unsigned int iy=0; iy < this->block_size; iy++) {
	rhs_holder[0]=E_y[0][iy];
	E_y[0][iy]-=C*dt/(2*dy) * (B_z[0][iy] - boundary_B_in[iy]);
	for(unsigned int ix=1; ix < this->block_size; ix++) {
	  rhs_holder[ix] = E_y[ix][iy];
	  E_y[ix][iy] -= C*dt/(2*dy) * (B_z[ix][iy] - B_z[ix-1][iy]);
	}
	for(unsigned int ix=0; ix < this->block_size-1; ix++) {
	  B_z[ix][iy] -= C*dt/(2*dy) * (rhs_holder[ix+1] - rhs_holder[ix]);
	}
	// if (y_line.rank() == 0 && iy == 0)
	//   for(unsigned int ix=0; ix < this->block_size; ix++)
	// 	std::cout << (ix == 0 ? "rhs_holder: " : "") << rhs_holder[ix] << (ix == block_size -1 ? "\n" : " ");
	  // std::cout << "BOUNDARY " << x_line.rank() << ": " << boundary_E_in[iy]
	  // 			<< " " << rhs_holder[this->block_size-1] << std::endl;
	B_z[this->block_size -1][iy] -= C*dt/(2*dy)*(boundary_E_in[iy] - rhs_holder[this->block_size-1]);
  }
}

void EMSquare2D::exchange_bdy_values(mpi::communicator comm, bdy_dir dir) {
  if(dir == BDY_X) {
	for(unsigned int iy=0; iy < this->block_size; iy++) {
	  boundary_E_out[iy] = E_y[0][iy];
	  boundary_B_out[iy] = B_z[this->block_size-1][iy];
	}
  } else {
	for(unsigned int ix=0; ix < this->block_size; ix++) {
	  boundary_E_out[ix] = E_x[ix][0];
	  boundary_B_out[ix] = B_z[ix][this->block_size-1];
	}
  }
  unsigned int p = comm.rank();
  unsigned int n_p = comm.size();
  if(p % 2 == 0) {
	if(p != n_p - 1)
	  comm.send(p+1,0,boundary_B_out, this->block_size);
	if(p != 0)
	  comm.recv(p-1,0, boundary_B_in, this->block_size);
	if(p != 0)
	  comm.send(p-1,0,boundary_E_out, this->block_size);
	if(p != n_p - 1)
	  comm.recv(p+1, 0, boundary_E_in, this->block_size);
  } else {
	if(p != 0)
	  comm.recv(p-1,0, boundary_B_in, this->block_size);
	if(p != n_p - 1)
	  comm.send(p+1,0,boundary_B_out, this->block_size);
	if(p != n_p - 1)
	  comm.recv(p+1,0, boundary_E_in, this->block_size);
	if(p != 0)
	  comm.send(p-1,0,boundary_E_out, this->block_size);
  }

  if(p == 0) {
	switch(dir) {
	case BDY_X:
	  for(unsigned int iy=0; iy < this->block_size; iy++) {
		this->boundary_B_in[iy] = this->B_z[0][iy];
	  }
	  break;
	case BDY_Y:
	  for(unsigned int ix=0; ix < this->block_size; ix++) {
		this->boundary_B_in[ix] = this->B_z[ix][0];
	  }
	}
  }
  if (p == n_p - 1) {
	std::fill_n(this->boundary_E_in, this->block_size, 0);
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

double EMSquare2DInitializer::x(unsigned int i) {
  return this->x_offset + (i + .5)*dx;
}

double EMSquare2DInitializer::y(unsigned int j) {
  return this->y_offset + (j + .5)*dy;
}

TE10Initializer::TE10Initializer(double dx, double dy, double L_x, double L_y,
								 unsigned int block_size, const mpi::communicator & comm)
: EMSquare2DInitializer(dx, dy, L_x, L_y, block_size, comm) {}

double TE10Initializer::E_x(unsigned int i, unsigned int j) {
  return 0;
}

double TE10Initializer::E_y(unsigned int i, unsigned int j) {
  return 0;
}

double TE10Initializer::B_z(unsigned int i, unsigned int j) {
  double x = this->x(i);

  return cos(M_PI*x/L_x);
}
