#include "AbstractReducedRHS.hpp"
#include <iostream>

RemoteReducedRHS::RemoteReducedRHS(TDCoupling* coupling, mpi::communicator& world, unsigned int block_size, unsigned int root) :
  AbstractReducedRHS(world, block_size) { 
  std::vector<mpi::request> requests(4);
  if(world.rank() != 0)
    requests[0] = world.isend(root, COUPLING_UPPER_1,
			      coupling->couplingCorrectionUpper[0]);
  if (world.rank() != world.size() - 1)
    requests[1] = world.isend(root, COUPLING_LOWER_M,
			      coupling->couplingCorrectionLower[this->blockSize-1]);
  if(world.rank() != 0 && world.rank() != world.size() -1) {
    requests[2] = world.isend(root, COUPLING_LOWER_1,
			      coupling->couplingCorrectionLower[0]);
    requests[3] = world.isend(root, COUPLING_UPPER_M,
			      coupling->couplingCorrectionUpper[this->blockSize-1]);
  }
  mpi::wait_all(requests.begin(), requests.end());
}
