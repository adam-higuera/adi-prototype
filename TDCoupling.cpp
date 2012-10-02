#include "AbstractReducedRHS.hpp"
#include "TDCoupling.hpp"

TDCoupling::TDCoupling(AbstractCouplingInitializer* init, LocalSolver* solver, unsigned int block_size) :
  couplingToLower(init->getCouplingToLower()),
  couplingFromLower(init->getCouplingFromLower()),
  couplingToUpper(init->getCouplingToUpper()),
  couplingFromUpper(init->getCouplingFromUpper()),
  couplingCorrectionLower(NULL),
  couplingCorrectionUpper(NULL),
  blockSize(block_size) {
  if(init->getBdyLoc() != LOWER_BDY && init->getBdyLoc() != DOUBLE_BDY) {
    couplingCorrectionLower = new double[blockSize];
    std::fill_n(couplingCorrectionLower, blockSize, 0);
    couplingCorrectionLower[blockSize-1] = couplingFromLower;
    solver->solve(couplingCorrectionLower);
  }
  if(init->getBdyLoc() != UPPER_BDY && init->getBdyLoc() != DOUBLE_BDY) {
    couplingCorrectionUpper = new double[blockSize];
    std::fill_n(couplingCorrectionUpper, blockSize, 0);
    couplingCorrectionUpper[0]= couplingFromUpper;
    solver->solve(couplingCorrectionUpper);
  }
}

void TDCoupling::applyCoupling(double* rhs, AbstractReducedRHS* red_rhs) {
  for(int i=0; i<this->blockSize; i++) {
    if(couplingCorrectionLower) {
      rhs[i] -= red_rhs->getLocalPart()[1]*this->couplingCorrectionLower[i];
    }
    if(couplingCorrectionUpper) {
      rhs[i] -= red_rhs->getLocalPart()[0]*this->couplingCorrectionUpper[i];
    }
  }  
}
