#include "LocalSolver.hpp"
#include <iostream>

LocalSolver::LocalSolver (AbstractMatrixInitializer* init, unsigned int block_size)
  : blockSize(block_size), diag(new double[blockSize]), upperDiag(new double[blockSize-1]),
    lowerDiag(new double[blockSize-1]), upperDiag2(new double[blockSize-1]), pivotPermutations(new int[blockSize]) {
  this->initializeMatrix(init);
}

void LocalSolver::initializeMatrix(AbstractMatrixInitializer* init) {
  int info;

  for (unsigned int i = 0; i < blockSize; i++) {
    this->diag [i] = (*init)(i, MAIN_DIAG);
    if (i < this->blockSize - 1) {
      this->lowerDiag[i] = (*init)(i, LOWER_DIAG);
      this->upperDiag [i] = (*init)(i, UPPER_DIAG);
    }
  }

  dgttrf_ (& this->blockSize, this->lowerDiag, this->diag,
		   this->upperDiag, this->upperDiag2, this->pivotPermutations, & info);
}

void LocalSolver::solve(double* rhs) {
  int one = 1;
  int info;

  dgttrs_("T", & this->blockSize, & one, this->lowerDiag, this->diag, this->upperDiag,
	  this->upperDiag2, this->pivotPermutations,
	  rhs,
	  & this->blockSize, & info);  
}
