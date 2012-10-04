#ifndef RHS_COMMUNICATOR_H
#define RHS_COMMUNICATOR_H

#include "common.hpp"
#include "AbstractReducedRHS.hpp"

class AbstractRHSCommunicator {
public:
  AbstractRHSCommunicator() {}
  AbstractRHSCommunicator(mpi::communicator& world) {}
  ~AbstractRHSCommunicator() {}

  virtual void doReducedSystems(std::vector<AbstractReducedRHS*> red_rhss)=0;
};

// class IndividualRHSCommunicator : public AbstractRHSCommunicator {
// public:
//   IndividualRHSCommunicator(mpi::communicator& world);
//   ~IndividualRHSCommunicator();

//   void doReducedSystems(std::vector<AbstractReducedRHS*> red_rhss);
// };

class CollectiveRHSCommunicator : public AbstractRHSCommunicator {
public:
  CollectiveRHSCommunicator(mpi::communicator& world, unsigned int block_size);
  ~CollectiveRHSCommunicator();

  void doReducedSystems(std::vector<AbstractReducedRHS*> red_rhss);

private:
};

#endif
