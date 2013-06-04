#ifndef ABSTRACT_REDUCED_RHS_H
#define ABSTRACT_REDUCED_RHS_H
#include "common.hpp"
#include "TDCoupling.hpp"

class AbstractReducedRHS {
public:
  AbstractReducedRHS(mpi::communicator& world, unsigned int block_size) : world(world), blockSize(block_size) {}
  ~AbstractReducedRHS() {}

  enum wich_matrix_element {COUPLING_UPPER_1, COUPLING_UPPER_M, COUPLING_LOWER_1, COUPLING_LOWER_M};

  virtual void copyValues(double* buf, unsigned int il, unsigned int n_l_ip) {}
  virtual void writeValues(double* buf, unsigned int il, unsigned int n_l_ip) {}
  virtual void solve() {}
  virtual double* getLocalPart() {return localPart;}

  virtual void yhallothar() {std::cout << "Y HALLO THAR" << std::endl;}

protected:
  unsigned int blockSize;
  double localPart[2];
  mpi::communicator& world;
};

class LocalReducedRHS : public AbstractReducedRHS {
public:
  LocalReducedRHS(TDCoupling* coupling, mpi::communicator& world, unsigned int block_size, unsigned int root);
  void copyValues(double* buf, unsigned int il, unsigned int n_l_ip);
  void writeValues(double* buf, unsigned int il, unsigned int n_l_ip);
  void solve();

  double* reducedDiag;
  double* reducedUpperDiag;
  double* reducedLowerDiag;
  double* reducedUpperDiag2; //For the factorization
  int* reducedPivotPermutations; //Ditto
  double* rhsStorage;

private:
  int one;
  int info;
  int reducedSize;
  unsigned int p;
};

class RemoteReducedRHS : public AbstractReducedRHS {
public:
  RemoteReducedRHS(TDCoupling* coupling, mpi::communicator& world, unsigned int block_size, unsigned int root);
};
#endif
