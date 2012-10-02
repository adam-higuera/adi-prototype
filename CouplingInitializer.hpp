#ifndef COUPLING_INITIALIZER_H
#define COUPLING_INITIALIZER_H
#include "common.hpp"
#include "MatrixInitializer.hpp"

class AbstractCouplingInitializer {
public:
  AbstractCouplingInitializer(AbstractMatrixInitializer* mat_init,
			      unsigned int block_size, mpi::communicator& world)
    : matInit(mat_init), blockSize(block_size), world(world) {}
  ~AbstractCouplingInitializer() {}

  virtual double getCouplingToLower() {
    return (*matInit)(blockSize-2, LOWER_DIAG);
  }
  virtual double getCouplingFromLower()=0;
  virtual double getCouplingToUpper() {
    return (*matInit)(0, UPPER_DIAG);
  };
  virtual double getCouplingFromUpper()=0;

  BoundaryLocation getBdyLoc() {
    if(world.size() == 1)
      return DOUBLE_BDY;
    if(world.rank() == 0)
      return UPPER_BDY;
    if(world.rank() == world.size() - 1)
      return LOWER_BDY;
    return NO_BDY;
  }

protected:
  AbstractMatrixInitializer* matInit;
  unsigned int blockSize;
  mpi::communicator& world;
};

class VacuumCouplingInitializer : public AbstractCouplingInitializer {
public:
  VacuumCouplingInitializer(AbstractMatrixInitializer* mat_init,
			    unsigned int block_size, mpi::communicator& world)
    : AbstractCouplingInitializer(mat_init, block_size, world) {}
  ~VacuumCouplingInitializer() {};

  // Vacuum Matrix is symmetric, no need to communicate
  double getCouplingFromLower() {
    return (*matInit)(blockSize-2, LOWER_DIAG);
  }
  double getCouplingFromUpper() {
    return (*matInit)(0, UPPER_DIAG);
  }
};

#endif
