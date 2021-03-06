#ifndef SIMULATION_H
#define SIMULATION_H

#include "RHSCollection.hpp"
#include <iostream>

class SimulationInitializer {
public:
  SimulationInitializer(double dx, double dy, double L_x, double L_y, unsigned int x_procs, unsigned int y_procs,
			unsigned int block_size, const mpi::communicator & comm);
  virtual double E_x(unsigned int i, unsigned int j) = 0;
  virtual double E_y(unsigned int i, unsigned int j) = 0;
  virtual double B_z(unsigned int i, unsigned int j) = 0;
  virtual AbstractRHSCollection* initCollection(std::vector<AbstractMatrixInitializer*> mat_inits,
						std::vector<AbstractCouplingInitializer*> c_inits,
						unsigned int block_size,
						mpi::communicator& comm)=0;
protected:
  double x_for_B(unsigned int i);
  double y_for_B(unsigned int j);
  double x_for_E(unsigned int i);
  double y_for_E(unsigned int j);
  
  unsigned int x_procs;
  unsigned int y_procs;

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
	     unsigned int n_cells, unsigned int x_procs, unsigned int y_procs, unsigned int n_steps,
	     unsigned int block_size,
	     std::string& dump_dir,
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

  double B_ul_corner;
  double B_ur_corner;
  double B_ll_corner;
  double B_lr_corner;

  double preFactorX;
  double preFactorY;

  double* BdyOutBufferTopRight;
  double* BdyOutBufferBotLeft;

  AbstractRHSCollection* xUpdateRHSs;
  AbstractRHSCollection* yUpdateRHSs;

  std::string dumpDir;

  void allocate_fields(SimulationInitializer* init);

  double getBdyB_z(mpi::communicator comm);
  void exchangeBdyValues(mpi::communicator comm, bdyDir dir);

  virtual void implicitUpdateM();
  virtual void explicitUpdateP();
  virtual void implicitUpdateP();
  virtual void explicitUpdateM();

  void implicitMSubstituteB();
  void implicitPSubstituteB();

  void yeeUpdate();

  void TimeStep();
  void printField(std::string msg);

  void dumpFields(std::string filename);
};

template<int m, int n, typename T>
class TEmnInitializer : public SimulationInitializer {
public:
  TEmnInitializer(double dx, double dy, double L_x, double L_y,
		  unsigned int x_procs, unsigned int y_procs,
		  unsigned int block_size, const mpi::communicator & comm)
    : SimulationInitializer(dx, dy, L_x, L_y, x_procs, y_procs, block_size, comm) {}
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
  virtual AbstractRHSCollection* initCollection(std::vector<AbstractMatrixInitializer*> mat_inits,
						std::vector<AbstractCouplingInitializer*> c_inits,
						unsigned int block_size,
						mpi::communicator& comm) {
    return new T(mat_inits, c_inits, block_size, comm);
  }
};

#endif
