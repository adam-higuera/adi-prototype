#ifndef LOCAL_SOLVER_H
#define LOCAL_SOLVER_H

#include "common.hpp"
#include "MatrixInitializer.hpp"

class LocalSolver {
public:
  LocalSolver(AbstractMatrixInitializer* init, unsigned int block_size);
  ~LocalSolver();

  void solve(double* rhs);
  void initializeMatrix(AbstractMatrixInitializer* init);

  int blockSize;
  double* diag;
  double* upperDiag;
  double* lowerDiag;
  double* upperDiag2; // Needed to store LU decomposition of local block
  int* pivotPermutations;
};
#endif
