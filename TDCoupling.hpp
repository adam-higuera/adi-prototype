#ifndef TD_COUPLING_H
#define TD_COUPLING_H

#include "common.hpp"
#include "CouplingInitializer.hpp"
#include "LocalSolver.hpp"

class AbstractReducedRHS;

class TDCoupling {
  friend class LocalReducedRHS;
  friend class RemoteReducedRHS;
public:
  TDCoupling(AbstractCouplingInitializer* init, LocalSolver* solver, unsigned int block_size);
  ~TDCoupling();

  void applyCoupling(double* RHS, AbstractReducedRHS* red_rhs);

  double couplingToLower;
  double couplingFromLower;
  double couplingToUpper;
  double couplingFromUpper;

private:
  unsigned int blockSize;
  
  double* couplingCorrectionUpper;
  double* couplingCorrectionLower;
};

#endif
