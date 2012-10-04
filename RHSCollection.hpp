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

class AbstractRHSCollection {
public:
  AbstractRHSCollection(std::vector<AbstractMatrixInitializer*> mat_inits,
						std::vector<AbstractCouplingInitializer*> coupling_inits,
						unsigned int block_size,
						mpi::communicator& world);
  void doLines(double** theLines);

  ReducedRHSFactory theFactory;
  unsigned int blockSize;
  mpi::communicator& world;
  unsigned int numLocalSolves;

  std::vector<TDCoupling*> couplings;
  std::vector<LocalSolver*> solvers;
  std::vector<AbstractReducedRHS*> redRHSs;

  double** rhsStorage;
};

class CollectiveRHSCollection : public AbstractRHSCollection {
public:
  RHSCollection(std::vector<AbstractMatrixInitializer*> mat_inits,
		std::vector<AbstractCouplingInitializer*> coupling_inits,
		unsigned int block_size,
		mpi::communicator& world);

  void doReducedSystems(std::vector<AbstractReducedRHS*> red_rhss);
  void dumpLine(unsigned int il, mpi::communicator& world);
private:
  double* sendbuf;
  double* recvbuf;
};

#endif
