#include "RHSCollection.hpp"

AbstractReducedRHS* ReducedRHSFactory::makeReducedRHS(unsigned int il, TDCoupling* coupling,
						      unsigned int block_size) {
  if(il % world.size() == world.rank())
    return new LocalReducedRHS(coupling, world, block_size, il % world.size());
  else
    return new RemoteReducedRHS(coupling, world, block_size, il % world.size());
}

AbstractRHSCollection::AbstractRHSCollection(std::vector<AbstractMatrixInitializer*> mat_inits,
					     std::vector<AbstractCouplingInitializer*> coupling_inits,
					     unsigned int block_size, mpi::communicator& world)
  :  world(world),
     blockSize(block_size),
     numLocalSolves(block_size / world.size() + (world.rank() < block_size % world.size())),
     couplings(block_size, NULL),
     solvers(block_size, NULL),
     redRHSs(block_size, NULL),
     rhsStorage(new double*[block_size]),
     theFactory(world)
{
  rhsStorage[0] = new double[blockSize*blockSize];
  for(unsigned int il=0; il < blockSize; il++) {
    if(il < blockSize-1)
      rhsStorage[il+1] = rhsStorage[il] + blockSize;
    solvers[il] = new LocalSolver(mat_inits[il], blockSize);
    couplings[il] = new TDCoupling(coupling_inits[il], solvers[il], blockSize);
    redRHSs[il] = theFactory.makeReducedRHS(il, couplings[il], blockSize);
  }  
}

CollectiveRHSCollection::CollectiveRHSCollection(std::vector<AbstractMatrixInitializer*> mat_inits,
						 std::vector<AbstractCouplingInitializer*> coupling_inits,
						 unsigned int block_size,
						 mpi::communicator& world)
  :  AbstractRHSCollection(mat_inits, coupling_inits, block_size, world),
	 sendbuf(new double[2*(block_size / world.size() + 1)]),
	 recvbuf(new double[2*world.size()*numLocalSolves])
 {}

void CollectiveRHSCollection::doLines(double** theLines) {
  for(unsigned int il=0; il < blockSize; il++) {
    solvers[il]->solve(theLines[il]);
    redRHSs[il]->getLocalPart()[0] = theLines[il][0];
    redRHSs[il]->getLocalPart()[1] = theLines[il][blockSize-1];
  }

  this->doReducedSystems(redRHSs);

  for(unsigned int il=0; il < blockSize; il++) {
    couplings[il]->applyCoupling(theLines[il], redRHSs[il]);
  }
}

void CollectiveRHSCollection::doReducedSystems(std::vector<AbstractReducedRHS*> red_rhss) {
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
    unsigned int n_l_ip = blockSize / world.size() + (ip < blockSize % world.size());
    mpi::scatter(world, recvbuf, sendbuf, 2*n_l_ip, ip);

    for(unsigned int il=0; il*world.size() + ip < blockSize; il++) {
      red_rhss[il*world.size() + ip]->getLocalPart()[0] = sendbuf[2*il];
      red_rhss[il*world.size() + ip]->getLocalPart()[1] = sendbuf[2*il+1];
    }
  }
}


void CollectiveRHSCollection::dumpLine(unsigned int il, mpi::communicator& world) {
  for(unsigned int ip=0; ip < world.size(); ip++) {
    if (world.rank() == ip)
      for(unsigned int i = 0; i < blockSize; i++)
	std::cout << rhsStorage[il][i] << " ";
    world.barrier();
  }
  std::cout << std::endl;
}
