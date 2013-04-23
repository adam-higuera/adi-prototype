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
  }
}


CollectiveRHSCollection::CollectiveRHSCollection(std::vector<AbstractMatrixInitializer*> mat_inits,
						 std::vector<AbstractCouplingInitializer*> coupling_inits,
						 unsigned int block_size,
						 mpi::communicator& world)
  : AbstractRHSCollection(mat_inits, coupling_inits, block_size, world),
    sendbuf(NULL),
    recvbuf(NULL) {
  if(world.size()==1)
    return;
  if(world.rank()==0) {
    recvbuf = new double[2*block_size*world.size()];
  }
  sendbuf = new double[2*block_size];
  for(unsigned int il=0; il < blockSize; il++) {
    if(world.rank() == 0)
      redRHSs[il] = new LocalReducedRHS(couplings[il], world, blockSize, 0);
    else
      redRHSs[il] = new RemoteReducedRHS(couplings[il], world, blockSize, 0);
  }
}

void CollectiveRHSCollection::doLines(double** theLines) {

#ifndef NO_LOCAL_SOLVES
  for(unsigned int il=0; il < blockSize; il++) {
    solvers[il]->solve(theLines[il]);
    if(world.size()!=1) {
      redRHSs[il]->getLocalPart()[0] = theLines[il][0];
      redRHSs[il]->getLocalPart()[1] = theLines[il][blockSize-1];
    }
  }
#endif
  if(world.size()==1)
    return;

#ifndef LOCAL_ONLY
  this->doReducedSystems(redRHSs);
#endif
#ifndef NO_APPLY_COUPLING
  for(unsigned int il=0; il < blockSize; il++) {
    couplings[il]->applyCoupling(theLines[il], redRHSs[il]);
  }
#endif
}

void CollectiveRHSCollection::doReducedSystems(std::vector<AbstractReducedRHS*>& red_rhss) {
#ifndef COMMUNICATION_ONLY
  for(int il=0; il < blockSize; il++) {
    std::memcpy(sendbuf + 2*il, red_rhss[il]->getLocalPart(), 2*sizeof(double));
  }
#endif
#ifndef NO_COLLECTIVES
  mpi::gather(world, sendbuf, 2*blockSize, recvbuf, 0);
#endif

#ifndef NO_REDUCED_SOLVE
  for(unsigned int il=0; il < blockSize; il++) {
    red_rhss[il]->copyValues(recvbuf, il, blockSize);
    red_rhss[il]->solve();
    red_rhss[il]->writeValues(recvbuf, il, blockSize);
  }
#endif

#ifndef NO_COLLECTIVES
  mpi::scatter(world, recvbuf, sendbuf, 2*blockSize, 0);
#endif
#ifndef COMMUNICATION_ONLY
  for(int il=0; il < blockSize; il++) {
    std::memcpy(red_rhss[il]->getLocalPart(), sendbuf + 2*il, 2*sizeof(double));
  }
#endif
}

ChunkedRHSCollection::ChunkedRHSCollection(std::vector<AbstractMatrixInitializer*> mat_inits,
						 std::vector<AbstractCouplingInitializer*> coupling_inits,
						 unsigned int block_size,
						 mpi::communicator& world)
  :  AbstractRHSCollection(mat_inits, coupling_inits, block_size, world),
	 sendbuf(new double[2*(block_size / world.size() + 1)]),
	 recvbuf(new double[2*world.size()*numLocalSolves]) {
  if(world.size()==1)
    return;
  for(unsigned int il=0; il < blockSize; il++)
    redRHSs[il] = (world.rank() == il % world.size()
		   ? (AbstractReducedRHS*) new LocalReducedRHS(couplings[il], world, blockSize, il % world.size())
		   : (AbstractReducedRHS*) new RemoteReducedRHS(couplings[il], world, blockSize, il % world.size()));
}

void ChunkedRHSCollection::doLines(double** theLines) {
#ifndef NO_LOCAL_SOLVES
  for(unsigned int il=0; il < blockSize; il++) {
    solvers[il]->solve(theLines[il]);
    if(world.size() != 1) {
      redRHSs[il]->getLocalPart()[0] = theLines[il][0];
      redRHSs[il]->getLocalPart()[1] = theLines[il][blockSize-1];
    }
  }
#endif

  if(world.size() == 1)
    return;

#ifndef LOCAL_ONLY
  this->doReducedSystems(redRHSs);
#endif

#ifndef COMMUNICATION_ONLY
  for(unsigned int il=0; il < blockSize; il++) {
    couplings[il]->applyCoupling(theLines[il], redRHSs[il]);
  }
#endif
}

void ChunkedRHSCollection::doReducedSystems(std::vector<AbstractReducedRHS*>& red_rhss) {
  unsigned int n_l_thisp = blockSize / world.size() + (world.rank() < blockSize % world.size());
  for(unsigned int ip=0; ip < world.size(); ip++) {
    // Number of reduced solves assigned to ip
    unsigned int n_l_ip = blockSize / world.size() + (ip < blockSize % world.size());
    // Get ready to send local part of reduced system assigned to processor ip
    if(n_l_ip == 0)
      break;

#ifndef COMMUNICATION_ONLY
    for(unsigned int il=0; il < n_l_ip; il++ ) {
      sendbuf[2*il] = red_rhss[il*world.size() + ip]->getLocalPart()[0];
      sendbuf[2*il+1] = red_rhss[il*world.size() + ip]->getLocalPart()[1];
    }
#endif    

#ifndef NO_COLLECTIVES
    mpi::gather(world, sendbuf, 2*n_l_ip, recvbuf, ip);
#endif
  }

#ifndef NO_REDUCED_SOLVE
  for(unsigned int il=0; il < blockSize; il++) {
    red_rhss[il]->copyValues(recvbuf, il, n_l_thisp);
    red_rhss[il]->solve();
    red_rhss[il]->writeValues(recvbuf, il, n_l_thisp);
  }
#endif

  for(unsigned int ip=0; ip < world.size(); ip++) {
    unsigned int n_l_ip = blockSize / world.size() + (ip < blockSize % world.size());
    if(n_l_ip == 0)
      break;

#ifndef NO_COLLECTIVES
    mpi::scatter(world, recvbuf, sendbuf, 2*n_l_ip, ip);
#endif

#ifndef COMMUNICATION_ONLY
    for(unsigned int il=0; il*world.size() + ip < blockSize; il++) {
      red_rhss[il*world.size() + ip]->getLocalPart()[0] = sendbuf[2*il];
      red_rhss[il*world.size() + ip]->getLocalPart()[1] = sendbuf[2*il+1];
    }
#endif
  }
}

void ChunkedRHSCollection::dumpLine(unsigned int il, mpi::communicator& world) {
  for(unsigned int ip=0; ip < world.size(); ip++) {
    if (world.rank() == ip)
      for(unsigned int i = 0; i < blockSize; i++)
	std::cout << rhsStorage[il][i] << " ";
    world.barrier();
  }
  std::cout << std::endl;
}

NonBlockingRHSCollection::NonBlockingRHSCollection(std::vector<AbstractMatrixInitializer*> mat_inits,
						   std::vector<AbstractCouplingInitializer*> coupling_inits,
						   unsigned int block_size,
						   mpi::communicator& world)
  : AbstractRHSCollection(mat_inits, coupling_inits, block_size, world),
    sendreqs(blockSize - numLocalSolves),
    scatter_recvreqs(blockSize - numLocalSolves),
    recvreqs(numLocalSolves, std::vector<mpi::request>(world.size()-1)),
    scatter_sendreqs(numLocalSolves, std::vector<mpi::request>(world.size()-1)),
    reduced_solve_buf(new double[numLocalSolves*2*(world.size()-1)]),
    recvbuf(new double[numLocalSolves*2*(world.size()-1)]),
    solves_done(numLocalSolves, false),
    corrections_done(blockSize, false) {}

void NonBlockingRHSCollection::doLines(double** theLines) {
  unsigned int p = world.rank();
  unsigned int n_p = world.size();
  unsigned int reduced_size = 2*(world.size() - 1);

  for(unsigned int il=0; il < blockSize; il++) {
    unsigned int proc_responsible = il % n_p;

    // Does another processor depend on this data?
    // If not, we'll compute it later
    if(proc_responsible != p) {
      solvers[il]->solve(theLines[il]);
      redRHSs[il]->getLocalPart()[0] = theLines[il][0];
      redRHSs[il]->getLocalPart()[1] = theLines[il][blockSize-1];
      
      // There should've been il/world.size() results we didn't have to send
      unsigned int send_num = il - il/n_p - (il % n_p > p);

      if(p == 0) {
	sendreqs[send_num] = world.isend(proc_responsible, il,
					 redRHSs[il]->getLocalPart()[1]);
      } else if (p == n_p - 1) {
	sendreqs[send_num] = world.isend(proc_responsible, il,
					 redRHSs[il]->getLocalPart()[0]);
      } else {
	sendreqs[send_num] = world.isend(proc_responsible, il,
					 redRHSs[il]->getLocalPart(), 2);
      }
    }
  }

  // Have now sent all data other processors are waiting for, so compute data
  // other processors don't depend on, and request data from other processors
  for(unsigned int il = p; il < blockSize; il += n_p) {
    solvers[il]->solve(theLines[il]);
    redRHSs[il]->getLocalPart()[0] = theLines[il][0];
    redRHSs[il]->getLocalPart()[1] = theLines[il][blockSize-1];

    unsigned int req_l_index = il / n_p;

    if(p != 0)
      reduced_solve_buf[reduced_size*req_l_index + 2*p - 1] = theLines[il][0];
    if(p != n_p - 1)
      reduced_solve_buf[reduced_size*req_l_index + 2*p] = theLines[il][blockSize-1];

    for(unsigned int ip = 0; ip < n_p; ip++) {
      unsigned int req_p_index = ip - (ip > p ? 1 : 0);
      if(ip != p) {
	if(ip == 0) {
	  recvreqs[req_l_index][req_p_index] = world.irecv(ip, il, reduced_solve_buf[0]);
	} else if (ip == n_p - 1) {
	  recvreqs[req_l_index][req_p_index] =
	    world.irecv(ip, il, reduced_solve_buf[reduced_size*(req_l_index+1) - 1]);
	} else {
	  recvreqs[req_l_index][req_p_index] =
	    world.irecv(ip, il, reduced_solve_buf + reduced_size*req_l_index + 2*ip - 1, 2);
	}
      }
    }
  }

  // Request results of reduced solves on other processors
  for(unsigned int il=0; il < blockSize; il++) {
    unsigned int recv_num = il - il/n_p - (il % n_p > p);
    if(il % n_p != p) {
      if(p == 0)
	scatter_recvreqs[recv_num] = world.irecv(il % n_p, il, redRHSs[il]->getLocalPart()[1]);
      else if(p == n_p - 1)
	scatter_recvreqs[recv_num] = world.irecv(il % n_p, il, redRHSs[il]->getLocalPart()[0]);
      else
	scatter_recvreqs[recv_num] = world.irecv(il % n_p, il, redRHSs[il]->getLocalPart(), 2);
    }
  }

  unsigned int red_solves_done = 0;
  unsigned int corrections_applied = 0;
  while (red_solves_done < numLocalSolves || corrections_applied < blockSize) {
    // Check to see whether all data required for reduced solves has been received
    // If so, perform reduced solve, send results to processors that need them,
    // and use results to apply correction locally.
    for(unsigned int il=0; il < numLocalSolves; il++) {
      if(solves_done[il] == false && mpi::test_all(recvreqs[il].begin(), recvreqs[il].end())) {
	unsigned int line_index = il*world.size() + world.rank();
	redRHSs[line_index]->copyValues(reduced_solve_buf + il*reduced_size - 1, 0, 1);
	redRHSs[line_index]->solve();
	redRHSs[line_index]->writeValues(reduced_solve_buf + il*reduced_size - 1, 0, 1);
	for(unsigned int ip=0; ip < world.size(); ip++) {
	  unsigned int req_p_index = ip - (ip > p ? 1 : 0);
	  if(ip != p) {
	    if(ip == 0)
	      scatter_sendreqs[il][req_p_index] = world.isend(ip, line_index,
						     reduced_solve_buf[il*reduced_size]);
	    else if(ip == n_p - 1)
	      scatter_sendreqs[il][req_p_index] = world.isend(ip, line_index,
						     reduced_solve_buf[(il+1)*reduced_size - 1]);
	    else
	      scatter_sendreqs[il][req_p_index] = world.isend(ip, line_index,
						     reduced_solve_buf + reduced_size*il + 2*ip - 1, 2);
	  } else {
	    if(p != 0)
	      redRHSs[line_index]->getLocalPart()[0] = reduced_solve_buf[reduced_size*il + 2*ip - 1];
	    if(p != n_p - 1)
	      redRHSs[line_index]->getLocalPart()[1] = reduced_solve_buf[reduced_size*il + 2*ip];
	    couplings[line_index]->applyCoupling(theLines[line_index], redRHSs[line_index]);
	    corrections_done[line_index] = true;
	  }
	}
	red_solves_done++;
	corrections_applied++;
	solves_done[il] = true;
      }
    }
    // Check to see whether other processors have sent data associated with their reduced solves,
    // if so, use data to apply local correction.
    for(unsigned int il=0; il < blockSize; il++) {
      unsigned int recv_num = il - il/n_p - (il % n_p > p);
      if(il % n_p != p) {
	if(corrections_done[il] == false && scatter_recvreqs[recv_num].test()) {
	  couplings[il]->applyCoupling(theLines[il], redRHSs[il]);
	  corrections_done[il] = true;
	  corrections_applied++;
	}
      }
    }
  }
  for(unsigned int il=0; il < numLocalSolves; il++) solves_done[il] = false;
  for(unsigned int il=0; il < blockSize; il++) corrections_done[il] = false;
}
