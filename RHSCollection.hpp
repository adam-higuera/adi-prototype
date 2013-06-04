#ifndef RHS_COLLECTION_H
#define RHS_COLLECTION_H
#include "common.hpp"
#include "LocalSolver.hpp"
#include "MatrixInitializer.hpp"
#include "CouplingInitializer.hpp"
#include "AbstractReducedRHS.hpp"

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
  virtual void doLines(double** theLines)=0;

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
  CollectiveRHSCollection(std::vector<AbstractMatrixInitializer*> mat_inits,
			  std::vector<AbstractCouplingInitializer*> coupling_inits,
			  unsigned int block_size,
			  mpi::communicator& world);

  void doReducedSystems(std::vector<AbstractReducedRHS*>& red_rhss);
  void doLines(double** theLines);
private:
  double* sendbuf;
  double* recvbuf;
};

class DelegatedRHSCollection : public AbstractRHSCollection {
public:
  DelegatedRHSCollection(std::vector<AbstractMatrixInitializer*> mat_inits,
			 std::vector<AbstractCouplingInitializer*> coupling_inits,
			 unsigned int block_size,
			 mpi::communicator& world);

  void doReducedSystems(std::vector<AbstractReducedRHS*>& red_rhss);
  void doLines(double** theLines);
private:
  double* sendbuf;
  double* recvbuf;
};


class ChunkedRHSCollection : public AbstractRHSCollection {
public:
  ChunkedRHSCollection(std::vector<AbstractMatrixInitializer*> mat_inits,
			  std::vector<AbstractCouplingInitializer*> coupling_inits,
			  unsigned int block_size,
			  mpi::communicator& world);

  void doReducedSystems(std::vector<AbstractReducedRHS*>& red_rhss);
  void dumpLine(unsigned int il, mpi::communicator& world);
  void doLines(double** theLines);
private:
  double* sendbuf;
  double* recvbuf;
};

class NonBlockingRHSCollection : public AbstractRHSCollection {
public:
  NonBlockingRHSCollection(std::vector<AbstractMatrixInitializer*> mat_inits,
			   std::vector<AbstractCouplingInitializer*> coupling_inits,
			   unsigned int block_size,
			   mpi::communicator& world);
  void doLines(double** theLines);
private:
  std::vector<mpi::request> sendreqs;
  std::vector<std::vector<mpi::request> > recvreqs;
  std::vector<mpi::request> scatter_recvreqs;
  std::vector<std::vector<mpi::request> > scatter_sendreqs;
  std::vector<bool> solves_done;
  std::vector<bool> corrections_done;
  double* reduced_solve_buf;
  double* recvbuf;
};

#endif
