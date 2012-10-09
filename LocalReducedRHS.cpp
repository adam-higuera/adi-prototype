#include "AbstractReducedRHS.hpp"

LocalReducedRHS::LocalReducedRHS(TDCoupling* coupling, mpi::communicator& world, unsigned int block_size, unsigned int root) :
  one(1),
  reducedSize(2*(world.size()-1)),
  p(world.rank()),
  AbstractReducedRHS(world, block_size),
  rhsStorage(new double[2*(world.size() - 1)]),
  reducedDiag(new double[2*(world.size()-1)]),
  reducedUpperDiag(new double[2*(world.size()-1)-1]),
  reducedLowerDiag(new double[2*(world.size()-1)-1]),
  reducedUpperDiag2(new double[2*(world.size()-1)-1]),
  reducedPivotPermutations(new int[2*(world.size()-1)-1]) {
  std::vector<mpi::request> requests(4*(world.size()-1));
  for(unsigned int ip=0; ip < world.size(); ip++) {
    unsigned int req_index = ip - (ip >= root);
    if (ip != root) {
      if (ip != 0) {
	requests[4*req_index] = world.irecv(ip, COUPLING_UPPER_1, this->reducedDiag[2*ip-1]);
      }
      if (ip != world.size() - 1) {
	requests[4*req_index+3] = world.irecv(ip, COUPLING_LOWER_M, this->reducedDiag[2*ip]);
      }
      if(ip != 0 && ip != world.size() - 1) {
	requests[4*req_index+1] = world.irecv(ip, COUPLING_UPPER_M, this->reducedLowerDiag[2*ip-1]);
	requests[4*req_index+2] = world.irecv(ip, COUPLING_LOWER_1, this->reducedUpperDiag[2*ip-1]);
      }
    } else {
      if (ip != 0)
	reducedDiag[2*ip-1] = coupling->couplingCorrectionUpper[0];
      if (ip != world.size()-1)
	reducedDiag[2*ip] = coupling->couplingCorrectionLower[blockSize-1];
      if(ip != 0 && ip != world.size()-1) {
	reducedLowerDiag[2*ip-1] = coupling->couplingCorrectionUpper[blockSize-1];
	reducedUpperDiag[2*ip-1] = coupling->couplingCorrectionLower[0];
      }
    }
  }

  for(unsigned int i=0; i < (world.size()-1); i++) {
    this->reducedUpperDiag[2*i] = this->reducedLowerDiag[2*i] = 1;
  }

  mpi::wait_all(requests.begin(), requests.end());

  dgttrf_ (& reducedSize, this->reducedLowerDiag, this->reducedDiag,
	   this->reducedUpperDiag, this->reducedUpperDiag2,
	   this->reducedPivotPermutations, & info);
}

void LocalReducedRHS::copyValues(double* buf, unsigned int il, unsigned int n_l_ip) {
  for(unsigned int ip=0; ip < world.size(); ip++) {
    if(ip != 0)
      rhsStorage[2*ip-1] = buf[2*n_l_ip*ip + 2*(il / world.size())];
    if(ip != world.size() - 1)
      rhsStorage[2*ip] = buf[2*n_l_ip*ip + 2*(il / world.size())+1];
  }
}

void LocalReducedRHS::writeValues(double* buf, unsigned int il, unsigned int n_l_ip) {
  for(unsigned int ip=0; ip < world.size(); ip++) {
    if(ip != 0)
      buf[2*n_l_ip*ip + 2*(il/world.size())] = rhsStorage[2*ip-1];
    if(ip != world.size() - 1)
      buf[2*n_l_ip*ip + 2*(il/world.size())+1] = rhsStorage[2*ip];
  }
}

void LocalReducedRHS::solve() {
  dgttrs_("T", & reducedSize, & one,
	  this->reducedLowerDiag, this->reducedDiag, this->reducedUpperDiag,
	  this->reducedUpperDiag2, this->reducedPivotPermutations,
	  this->rhsStorage,
	  & reducedSize, & info);
}
