#ifndef MATRIX_BLOCK_H
#define MATRIX_BLOCK_H

#include<boost/mpi.hpp>

extern "C" {
  int dgttrf_(int*, double*, double*, double*, double*, int*, int*);
  int dgttrs_(char*, int*, int*, double*, double*, double*, double*, int*, double*, int*, int*);
}

namespace mpi = boost::mpi;
enum which_diagonal {MAIN_DIAG, UPPER_DIAG, LOWER_DIAG};
enum wich_matrix_element {COUPLING_UPPER_1, COUPLING_UPPER_M, COUPLING_LOWER_1, COUPLING_LOWER_M};

class MatrixBlock {
public:
  template<typename Initializer>
  MatrixBlock (mpi::communicator & world, unsigned int problem_size, Initializer init);
  ~MatrixBlock();

  void solve(double* rhs);

  void printDiag();

private:

  void allocate_storage ();

  template<typename Initializer>
  void initialize_matrix(Initializer init);

  unsigned int global_i_from_local (int i) {
	return this->index_offset + i;
  }
  int block_size;
  unsigned int index_offset;

  void exchange_couplings();
  double coupling_to_lower, coupling_to_upper;
  double coupling_from_lower, coupling_from_upper;

  void compute_coupling_corrections();

  void gather_reduced_matrix();
  
  double* coupling_correction_lower;
  double* coupling_correction_upper;

  double* orig_diag;
  
  double* diag;
  double* upper_diag;
  double* lower_diag;
  double* U_upper_diag2; // Needed to store LU decomposition of local block
  int* pivot_permutations;

  double* reduced_diag;
  double* reduced_lower_diag;
  double* reduced_upper_diag;
  double* reduced_U_upper_diag2;
  int* reduced_pivot_permutations;

  double* reduced_rhs;

  const mpi::communicator & world;
};

class ToeplitzMatrixInitializer {
public:
  ToeplitzMatrixInitializer(double diag, double off_diag) : diag (diag), off_diag (off_diag) {}
  double operator ()(unsigned int i, which_diagonal w) const {
	switch (w) {
	case MAIN_DIAG:
	  return this->diag;
	case LOWER_DIAG:
	case UPPER_DIAG:
	  return this->off_diag;
	}
  }
private:
  double diag, off_diag;
};

class EMToeplitzMatrixInitializer {
public:
  EMToeplitzMatrixInitializer(double diag, double off_diag, unsigned int problem_size)
    : diag (diag), off_diag (off_diag), problem_size(problem_size) {}
  double operator ()(unsigned int i, which_diagonal w) const {
	switch (w) {
	case MAIN_DIAG:
	  return (i == 0 || i == problem_size - 1) ? this->diag + this->off_diag : this->diag;
	case LOWER_DIAG:
	case UPPER_DIAG:
	  return this->off_diag;
	}
  }
private:
  double diag, off_diag;
  unsigned int problem_size;
};

template<typename Initializer>
MatrixBlock::MatrixBlock (mpi::communicator & world, unsigned int problem_size, Initializer init)
  : diag(NULL), upper_diag(NULL), lower_diag(NULL), world(world) {
  this->block_size = problem_size / world.size () + (world.rank () < problem_size % world.size ());
  this->index_offset = world.rank() * (problem_size / world.size ())
	+ std::min(world.rank(), int(problem_size % world.size ()));

  this->allocate_storage();
  this->initialize_matrix(init);
  this->exchange_couplings();
  this->compute_coupling_corrections();

  this->gather_reduced_matrix();
}

template<typename Initializer>
void MatrixBlock::initialize_matrix(Initializer init) {
  int info;
  for (unsigned int i = 0; i < this->block_size; i++) {
    this->diag [i] = init (global_i_from_local(i), MAIN_DIAG);
    if (i < this->block_size - 1) {
      this->lower_diag [i] = init (this->global_i_from_local(i), LOWER_DIAG);
      this->upper_diag [i] = init (this->global_i_from_local(i), UPPER_DIAG);
    }
  }
  this->coupling_to_upper = init (this->global_i_from_local (-1), UPPER_DIAG);
  this->coupling_to_lower = init (this->global_i_from_local (this->block_size-1), LOWER_DIAG);

  std::copy(diag, diag + block_size, orig_diag);

  dgttrf_ (& this->block_size, this->lower_diag, this->diag,
		   this->upper_diag, this->U_upper_diag2, this->pivot_permutations, & info);
}


#endif
