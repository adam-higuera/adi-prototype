#include "RHSCommunicator.hpp"
#include <iostream>

CollectiveRHSCommunicator::CollectiveRHSCommunicator(mpi::communicator& world, unsigned int block_size) :
  world(world),
  blockSize(block_size),
  numLocalSolves(block_size / world.size() + (world.rank() < block_size % world.size())),
  sendbuf(new double[2*(block_size / world.size() + 1)]),
  recvbuf(new double[2*world.size()*numLocalSolves]) {}

void CollectiveRHSCommunicator::doReducedSystems(std::vector<AbstractReducedRHS*> red_rhss) {
  unsigned int n_l_thisp = blockSize / world.size() + (world.rank() < blockSize % world.size());
  for(unsigned int ip=0; ip < world.size(); ip++) {
    // Number of reduced solves assigned to ip
    unsigned int n_l_ip = blockSize / world.size() + (ip < blockSize % world.size());
    // Get ready to send local part of reduced system assigned to processor ip
    for(unsigned int il=0; il < n_l_ip; il++ ) {
      sendbuf[2*il] = red_rhss[il*world.size() + ip]->getLocalPart()[0];
      sendbuf[2*il+1] = red_rhss[il*world.size() + ip]->getLocalPart()[1];
    }

    mpi::gather(world, sendbuf, 2*n_l_ip, recvbuf, ip);
  }

  for(unsigned int il=0; il < blockSize; il++) {
    red_rhss[il]->copyValues(recvbuf, il, n_l_thisp);
    red_rhss[il]->solve();
    red_rhss[il]->writeValues(recvbuf, il, n_l_thisp);
  }

  for(unsigned int ip=0; ip < world.size(); ip++) {
    mpi::scatter(world, recvbuf, sendbuf, 2*numLocalSolves, ip);

    for(unsigned int il=0; il*world.size() + ip < blockSize; il++) {
      red_rhss[il*world.size() + ip]->getLocalPart()[0] = sendbuf[2*il];
      red_rhss[il*world.size() + ip]->getLocalPart()[1] = sendbuf[2*il+1];
    }
  }
}
