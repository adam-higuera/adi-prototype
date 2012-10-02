#ifndef RHS_COLLECTION_H
#define RHS_COLLECTION_H
#include "common.hpp"
#include "LocalSolver.hpp"
#include "MatrixInitializer.hpp"
#include "CouplingInitializer.hpp"
#include "RHSCommunicator.hpp"

class ReducedRHSFactory {
public:
  ReducedRHSFactory(mpi::communicator& world) : world(world) {}
  AbstractReducedRHS* makeReducedRHS(unsigned int il, TDCoupling* coupling, unsigned int block_size);
private:
  mpi::communicator& world;
};

class RHSCollection {
public:
  RHSCollection(std::vector<AbstractMatrixInitializer*> mat_inits,
		std::vector<AbstractCouplingInitializer*> coupling_inits,
		AbstractRHSCommunicator* the_comm,
		unsigned int block_size,
		mpi::communicator& world);
  void doLines(double** theLines);
  void dumpLine(unsigned int il, mpi::communicator& world);
private:
  ReducedRHSFactory theFactory;
  
  std::vector<TDCoupling*> couplings;
  std::vector<LocalSolver*> solvers;
  std::vector<AbstractReducedRHS*> redRHSs;
  AbstractRHSCommunicator* theComm;
  unsigned int blockSize;

  double** rhsStorage;
};

#endif
