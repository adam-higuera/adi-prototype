#include <iostream>
#include <algorithm>
#include "MatrixBlock.hpp"

MatrixBlock::~MatrixBlock () {
  delete[] this->diag;
  delete[] this->upper_diag;
  delete[] this->lower_diag;
}

void MatrixBlock::compute_coupling_corrections() {
  int one = 1;
  int info;
  if(world.rank() != world.size()-1) {
	std::fill_n(this->coupling_correction_lower, this->block_size, 0);
	this->coupling_correction_lower[block_size-1] = this->coupling_from_lower;
	dgttrs_("T", & this->block_size, & one, this->lower_diag, this->diag, this->upper_diag,
		   this->U_upper_diag2, this->pivot_permutations,
		   this->coupling_correction_lower,
			& this->block_size, & info);
  }
  if(world.rank() != 0) {
	std::fill_n(this->coupling_correction_upper, this->block_size, 0);
	this->coupling_correction_upper[0]=this->coupling_from_upper;
	dgttrs_("T", & this->block_size, & one, this->lower_diag, this->diag, this->upper_diag,
		   this->U_upper_diag2, this->pivot_permutations,
		   this->coupling_correction_upper,
			& this->block_size, & info);
  }
}

void MatrixBlock::gather_reduced_matrix() {
  int info;
  if(world.rank() != 0) {
	std::vector<mpi::request> requests(4);
	requests[0] = world.isend(0, COUPLING_UPPER_1, this->coupling_correction_upper[0]);
	if (world.rank() != world.size() - 1) {
	  requests[1] = world.isend(0, COUPLING_UPPER_M, this->coupling_correction_upper[this->block_size-1]);
	  requests[2] = world.isend(0, COUPLING_LOWER_1, this->coupling_correction_lower[0]);
	  requests[3] = world.isend(0, COUPLING_LOWER_M, this->coupling_correction_lower[this->block_size-1]);
	}
	mpi::wait_all(requests.begin(), requests.end());
  } else {
	std::vector<mpi::request> requests(4*(world.size()-2)+1);
	for(int i=1; i < world.size(); i++) {
	  requests[4*(i-1)] = world.irecv(i, COUPLING_UPPER_1, this->reduced_diag[2*i-1]);
	  if(i != world.size() - 1) {
		requests[4*(i-1)+1] = world.irecv(i, COUPLING_UPPER_M, this->reduced_lower_diag[2*i-1]);
		requests[4*(i-1)+2] = world.irecv(i, COUPLING_LOWER_1, this->reduced_upper_diag[2*i-1]);
		requests[4*(i-1)+3] = world.irecv(i, COUPLING_LOWER_M, this->reduced_diag[2*i]);
	  }
	  this->reduced_diag[0] = this->coupling_correction_lower[this->block_size-1];
	  for(int i=0; i < (world.size()-1); i++) {
		this->reduced_upper_diag[2*i] = this->reduced_lower_diag[2*i] = 1;
	  }
	}
	mpi::wait_all(requests.begin(), requests.end());
  }

  if(world.rank() == 0) {
	int reduced_size = 2*(world.size()-1);
	dgttrf_ (& reduced_size, this->reduced_lower_diag, this->reduced_diag,
			this->reduced_upper_diag, this->reduced_U_upper_diag2,
			 this->reduced_pivot_permutations, & info);	
  }
}

void MatrixBlock::allocate_storage() {
  this->diag = new double [this->block_size];
  this->orig_diag = new double[this->block_size];
  this->upper_diag = new double [this->block_size-1];
  this->lower_diag = new double [this->block_size-1];
  this->U_upper_diag2 = new double [this->block_size-1];
  this->pivot_permutations = new int [this->block_size-1];
  if(world.rank() != world.size() - 1) {
	this->coupling_correction_lower = new double[this->block_size];
  }
  if(world.rank() != 0) {
	this->coupling_correction_upper = new double[this->block_size];
  }

  if(world.rank() == 0) {
	this->reduced_diag = new double[2*(world.size()-1)];
	this->reduced_lower_diag = new double[2*(world.size()-1) - 1];
	this->reduced_upper_diag = new double[2*(world.size()-1) - 1];
	this->reduced_U_upper_diag2 = new double[2*(world.size()-1) - 1];
	this->reduced_pivot_permutations = new int[2*(world.size()-1)];
	// Has a zero of padding on either side so that we can use mpi::gather
	this->reduced_rhs = new double[2*world.size()];
  }
}

void MatrixBlock::exchange_couplings() {
  // Handle coupling between pairs of blocks - if even-numbered, exchange with block above, otherwise, with block below
  if (world.rank () % 2 == 0) {
	if (world.rank () != 0) {
	  world.send (world.rank ()-1, 0, this->coupling_to_upper);
	  world.recv (world.rank ()-1, 0, this->coupling_from_upper);
	}
	if(world.rank() != world.size() - 1) {
	  world.send(world.rank()+1, 0, this->coupling_to_lower);
	  world.recv(world.rank()+1, 0, this->coupling_from_lower);
	}
  } else {
	if(world.rank() != world.size()-1) {
	  world.recv(world.rank()+1, 0, this->coupling_from_lower);
	  world.send(world.rank()+1, 0, this->coupling_to_lower);
	}
	if(world.rank() != 0) {
	  world.recv(world.rank() -1, 0, this->coupling_from_upper);
	  world.send(world.rank()-1, 0, this->coupling_to_upper);
	}
  }
}

void MatrixBlock::solve(double* rhs) {
  int one = 1;
  int info;

  dgttrs_("T", & this->block_size, & one, this->lower_diag, this->diag, this->upper_diag,
		 this->U_upper_diag2, this->pivot_permutations,
		 rhs,
		  & this->block_size, & info);
  double values[2] = {0, 0};

  // Gather reduced system on process 0, solve it, scatter it
  if(world.rank() == 0) {
	int reduced_size = 2*(world.size()-1);
	values[1] = rhs[this->block_size-1];
	  
	mpi::gather(world, values, 2, this->reduced_rhs, 0);

	dgttrs_("T", & reduced_size, & one, this->reduced_lower_diag, this->reduced_diag, this->reduced_upper_diag,
		   this->reduced_U_upper_diag2, this->reduced_pivot_permutations,
		   this->reduced_rhs + 1, //There's a leading zero of padding
			& reduced_size, & info);

	mpi::scatter(world, this->reduced_rhs, values, 2, 0);
	  } else if (world.rank() == world.size()-1) {
	values[0] = rhs[0];

	mpi::gather(world, values, 2, 0);
	mpi::scatter(world, values, 2, 0);
  } else {
	values[0] = rhs[0];
	values[1] = rhs[this->block_size-1];

	mpi::gather(world, values, 2, 0);
	mpi::scatter(world, values, 2, 0);
  }

  for(int i=0; i<this->block_size; i++) {
	if(world.rank() != world.size()-1) {
	  rhs[i] -= values[1]*this->coupling_correction_lower[i];
	}
	if(world.rank() != 0) {
	  rhs[i] -= values[0]*this->coupling_correction_upper[i];
	}
  }
}

void MatrixBlock::printDiag() {
  int dummy;
  if(world.rank() == 0)
    std::cout << "DIAG: ";
  if(world.rank() != 0)
    world.recv(world.rank() - 1, 0, dummy);
  for(unsigned int ix=0; ix < block_size; ix++) {
    std::cout << orig_diag[ix] << " ";
  }
  std::cout.flush();
  if (world.rank() == world.size() - 1)
    std::cout << std::endl;
  else
    world.send(world.rank() + 1, 0, dummy);
}
