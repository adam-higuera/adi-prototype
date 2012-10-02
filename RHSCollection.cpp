#include "RHSCollection.hpp"

AbstractReducedRHS* ReducedRHSFactory::makeReducedRHS(unsigned int il, TDCoupling* coupling,
						      unsigned int block_size) {
  if(il % world.size() == world.rank())
    return new LocalReducedRHS(coupling, world, block_size, il % world.size());
  else
    return new RemoteReducedRHS(coupling, world, block_size, il % world.size());
}

RHSCollection::RHSCollection(std::vector<AbstractMatrixInitializer*> mat_inits,
			     std::vector<AbstractCouplingInitializer*> coupling_inits,
			     AbstractRHSCommunicator* the_comm,
			     unsigned int block_size,
			     mpi::communicator& world)
  : blockSize(block_size),
    couplings(block_size, NULL),
    solvers(block_size, NULL),
    redRHSs(block_size, NULL),
    rhsStorage(new double*[block_size]),
    theComm(the_comm),
    theFactory(world) {
  rhsStorage[0] = new double[blockSize*blockSize];
  for(unsigned int il=0; il < blockSize; il++) {
    if(il < blockSize-1)
      rhsStorage[il+1] = rhsStorage[il] + blockSize;
    solvers[il] = new LocalSolver(mat_inits[il], blockSize);
    couplings[il] = new TDCoupling(coupling_inits[il], solvers[il], blockSize);
    redRHSs[il] = theFactory.makeReducedRHS(il, couplings[il], blockSize);
  }
}

void RHSCollection::doLines(double** theLines) {
  for(unsigned int il=0; il < blockSize; il++) {
    solvers[il]->solve(theLines[il]);
    redRHSs[il]->getLocalPart()[0] = theLines[il][0];
    redRHSs[il]->getLocalPart()[1] = theLines[il][blockSize-1];
  }

  theComm->doReducedSystems(redRHSs);

  for(unsigned int il=0; il < blockSize; il++) {
    // for(int ip = 0; ip < blockSize; ip++)
    //   std::cout << theLines[il][ip] << " ";
    // std::cout << std::endl;
    couplings[il]->applyCoupling(theLines[il], redRHSs[il]);
    // for(int ip = 0; ip < blockSize; ip++)
    //   std::cout << theLines[il][ip] << " ";
    // std::cout << std::endl;
  }
}

void RHSCollection::dumpLine(unsigned int il, mpi::communicator& world) {
  for(unsigned int ip=0; ip < world.size(); ip++) {
    if (world.rank() == ip)
      for(unsigned int i = 0; i < blockSize; i++)
	std::cout << rhsStorage[il][i] << " ";
    world.barrier();
  }
  std::cout << std::endl;
}
