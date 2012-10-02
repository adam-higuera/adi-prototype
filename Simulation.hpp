#ifndef SIMULATION_H
#define SIMULATION_H

#include "RHSCommunicator.hpp"
#include "RHSCollection.hpp"
#include <iostream>

class SimulationInitializer {
public:
  SimulationInitializer(double dx, double dy, double L_x, double L_y,
			unsigned int block_size, const mpi::communicator & comm);
  virtual double E_x(unsigned int i, unsigned int j) = 0;
  virtual double E_y(unsigned int i, unsigned int j) = 0;
  virtual double B_z(unsigned int i, unsigned int j) = 0;
protected:
  double x_for_B(unsigned int i);
  double y_for_B(unsigned int j);
  double x_for_E(unsigned int i);
  double y_for_E(unsigned int j);
  
  double dx;
  double dy;
  double x_offset;
  double y_offset;
  double L_x;
  double L_y;
};

class Simulation {
public:
  Simulation(double L_x, double L_y, double T,
	     unsigned int n_cells, unsigned int n_steps,
	     SimulationInitializer* init, mpi::communicator & world);

  void simulate(bool dump=true, unsigned int dump_periodicity=9, unsigned int total_dumps=100);

protected:
  mpi::communicator & world, xLine, yLine;

  unsigned int blockSize, nSteps;

  double** E_x;
  double** E_y;
  double** B_z;

  double** rhss;

  double dx;
  double dy;
  double dt;

  double* BLeftBdy;
  double* BRightBdy;
  double* BBotBdy;
  double* BTopBdy;

  double* BdyOutBufferTopRight;
  double* BdyOutBufferBotLeft;

  RHSCollection* xUpdateRHSs;
  RHSCollection* yUpdateRHSs;

  void allocate_fields(SimulationInitializer* init);

  double getBdyB_z(mpi::communicator comm);
  void exchangeBdyValues(mpi::communicator comm, bdyDir dir);

  virtual void implicitUpdateM();
  virtual void explicitUpdateP();
  virtual void implicitUpdateP();
  virtual void explicitUpdateM();

  void implicitMSubstituteB();
  void implicitPSubstituteB();

  void TimeStep();
  void printField(std::string msg);

  void dumpFields(std::string filename);
};

template<int m, int n>
class TEmnInitializer : public SimulationInitializer {
public:
  TEmnInitializer(double dx, double dy, double L_x, double L_y,
		  unsigned int block_size, const mpi::communicator & comm)
    : SimulationInitializer(dx, dy, L_x, L_y, block_size, comm) {}
  virtual double E_x(unsigned int i, unsigned int j) {
    return 0;
  }
  virtual double E_y(unsigned int i, unsigned int j) {
    return 0;
  }
  virtual double B_z(unsigned int i, unsigned int j) {
    double x = this->x_for_B(i);
    double y = this->y_for_B(j);

    return cos(m*M_PI*x/L_x)*cos(n*M_PI*y/L_y);
  }
};

#endif
