#ifndef SIMULATION_3D_H
#define SIMULATION_3D_H

#include "RHSCollection.hpp"
#include "hdf5.h"
#include <iostream>

class Simulation3DInitializer {
public:
  Simulation3DInitializer(double dx, double dy, double dz,
			  double L_x, double L_y, double L_z,
			  unsigned int block_size);
  virtual void populate(double* E, double* B, unsigned int ix, unsigned int iy, unsigned int iz) = 0;
  virtual void populateGuard(double* guard_ptr_E, double* guard_ptr_B, int ix, int iy, int iz) = 0;
  virtual AbstractRHSCollection* initCollection(std::vector<AbstractMatrixInitializer*> mat_inits,
						std::vector<AbstractCouplingInitializer*> c_inits,
						unsigned int block_size,
						mpi::communicator& comm)=0;
  void setOffsets(const mpi::communicator & xLine, const mpi::communicator & yLine, const mpi::communicator & zLine);

  unsigned int blockSize;
  
  double dx;
  double dy;
  double dz;
  double x_offset;
  double y_offset;
  double z_offset;
  double L_x;
  double L_y;
  double L_z;
};


enum commDir {UPWARDS, DOWNWARDS};

class Simulation3D {
public:

  Simulation3D(double L_x, double L_y, double L_z,
	       double T,
	       unsigned int n_cells, unsigned int n_steps,
	       unsigned int procs_x, unsigned int procs_y, unsigned int procs_z,
	       unsigned int block_size,
	       std::string& dump_dir,
	       Simulation3DInitializer* init, mpi::communicator & world);
  ~Simulation3D();
  
  double dx;
  double dy;
  double dz;
  double dt;

  unsigned int procsX;
  unsigned int procsY;
  unsigned int procsZ;

  unsigned int blockSize;
  unsigned int nSteps;
  unsigned int currentStep;

  mpi::communicator xLine, yLine, zLine;
  mpi::communicator& world;

  double*** guardB;
  double*** guardE;

  double* E;
  double* B;
  double* tmp_field;

  double* rhsx;
  double** rhs_ptrs_x;
  double* rhsy;
  double** rhs_ptrs_y;
  double* rhsz;
  double** rhs_ptrs_z;
  double* guardSendbuf;
  
  double preFactorX;
  double preFactorY;
  double preFactorZ;

  AbstractRHSCollection* xUpdateRHSs;
  AbstractRHSCollection* yUpdateRHSs;
  AbstractRHSCollection* zUpdateRHSs;

  std::string dumpDir;

  void simulate(bool dump, unsigned int dump_periodicity, unsigned int total_dumps);
  double*** allocateGuardStorage();
  void initFields(Simulation3DInitializer* init);

  void fillSendbufX(unsigned int ix, double* F);
  void fillSendbufY(unsigned int iy, double* F);
  void fillSendbufZ(unsigned int iz, double* F);

  void getGuardF(double* F, double*** guardF, commDir dir);

  virtual void implicitUpdateM();
  virtual void explicitUpdateP();
  virtual void implicitUpdateP();
  virtual void explicitUpdateM();

  void populateRHSM();
  void populateRHSP();
  void writeRHSM();
  void writeRHSP();

  void implicitMSubstituteB();
  void implicitPSubstituteB();

  void yeeUpdate();

  void timeStep();
  // void printField(std::string msg);

  void dumpFields(std::string filename);
  void dumpTimings(unsigned long* timings, hsize_t total_timings, unsigned int steps_per_timing);
  void dumpCenterFields(double* fields, hsize_t nSteps);
};

template<int m, int n, int l, typename T>
class TEmnlInitializer : public Simulation3DInitializer {
public:
  TEmnlInitializer(double dx, double dy, double dz,
		   double L_x, double L_y, double L_z,
		   unsigned int block_size)
    : Simulation3DInitializer(dx, dy, dz, L_x, L_y, L_z, block_size),
      k_x(m*M_PI/L_x), k_y(n*M_PI/L_y), k_z(l*M_PI/L_z) {}

  virtual void populate(double* E, double* B,
			  unsigned int ix, unsigned int iy, unsigned int iz) {
    double x_E = x_offset + ix*dx; double x_B = x_offset + (ix + 0.5)*dx;
    double y_E = y_offset + iy*dx; double y_B = y_offset + (iy + 0.5)*dy;
    double z_E = z_offset + iz*dz; double z_B = z_offset + (iz + 0.5)*dz;

    unsigned int field_offset = ((blockSize*ix + iy)*blockSize + iz)*3;

    double* E_ptr = E + field_offset; double* B_ptr = B + field_offset;
    E_ptr[0] = 0; E_ptr[1] = 0; E_ptr[2] = 0;
    B_ptr[0] = (-k_x/k_y)*k_z*sin(k_x*x_B)*cos(k_y*y_B)*cos(k_z*z_B);
    B_ptr[1] = -k_z*cos(k_x*x_B)*sin(k_y*y_B)*cos(k_z*z_B);
    B_ptr[2] = (k_x/k_y*k_x + k_y)*cos(k_x*x_B)*cos(k_y*y_B)*sin(k_z*z_B);
  }

  virtual void populateGuard(double* guard_ptr_E, double* guard_ptr_B,
			     int ix, int iy, int iz) {
    double x_E = x_offset + ix*dx; double x_B = x_offset + (ix + 0.5)*dx;
    double y_E = y_offset + iy*dx; double y_B = y_offset + (iy + 0.5)*dy;
    double z_E = z_offset + iz*dz; double z_B = z_offset + (iz + 0.5)*dz;

    guard_ptr_E[0] = 0; guard_ptr_E[1] = 0; guard_ptr_E[2] = 0;
    guard_ptr_B[0] = (-k_x/k_y)*k_z*sin(k_x*x_B)*cos(k_y*y_B)*cos(k_z*z_B);
    guard_ptr_B[1] = -k_z*cos(k_x*x_B)*sin(k_y*y_B)*cos(k_z*z_B);
    guard_ptr_B[2] = (k_x/k_y*k_x + k_y)*cos(k_x*x_B)*cos(k_y*y_B)*sin(k_z*z_B);
  }

  virtual AbstractRHSCollection* initCollection(std::vector<AbstractMatrixInitializer*> mat_inits,
						std::vector<AbstractCouplingInitializer*> c_inits,
						unsigned int block_size,
						mpi::communicator& comm) {
    return new T(mat_inits, c_inits, block_size, comm);
  }

  double k_x;
  double k_y;
  double k_z;
};

template<int m, int n, int l, typename T>
class ShiftedTEmnlInitializer : public Simulation3DInitializer {
public:
  ShiftedTEmnlInitializer(double dx, double dy, double dz,
			  double L_x, double L_y, double L_z,
			  unsigned int block_size)
    : Simulation3DInitializer(dx, dy, dz, L_x, L_y, L_z, block_size),
      k_x(m*M_PI/L_x), k_y(n*M_PI/L_y), k_z(l*M_PI/L_z) {}

  virtual void populate(double* E, double* B,
			unsigned int ix, unsigned int iy, unsigned int iz) {
    double x_E = x_offset + ix*dx; double x_B = x_offset + (ix + 0.5)*dx;
    double y_E = y_offset + iy*dx; double y_B = y_offset + (iy + 0.5)*dy;
    double z_E = z_offset + iz*dz; double z_B = z_offset + (iz + 0.5)*dz;

    unsigned int field_offset = ((blockSize*ix + iy)*blockSize + iz)*3;

    double* E_ptr = E + field_offset; double* B_ptr = B + field_offset;
    B_ptr[0] = 0; B_ptr[1] = 0; B_ptr[2] = 0;
    E_ptr[0] = cos(k_x*(x_E + 0.5*dx))*sin(k_y*y_E)*sin(k_z*z_E);
    E_ptr[1] = -k_x/k_y*sin(k_x*x_E)*cos(k_y*(y_E + 0.5*dy))*sin(k_z*z_E);
    E_ptr[2] = 0;
  }

  virtual void populateGuard(double* guard_ptr_E, double* guard_ptr_B,
			     int ix, int iy, int iz) {
    double x_E = x_offset + ix*dx; double x_B = x_offset + (ix + 0.5)*dx;
    double y_E = y_offset + iy*dx; double y_B = y_offset + (iy + 0.5)*dy;
    double z_E = z_offset + iz*dz; double z_B = z_offset + (iz + 0.5)*dz;

    guard_ptr_B[0] = 0; guard_ptr_B[1] = 0; guard_ptr_B[2] = 0;
    guard_ptr_E[0] = cos(k_x*(x_E + 0.5*dx))*sin(k_y*y_E)*sin(k_z*z_E);
    guard_ptr_E[1] = -k_x/k_y*sin(k_x*x_E)*cos(k_y*(y_E + 0.5*dy))*sin(k_z*z_E);
    guard_ptr_E[2] = 0;
  }

  virtual AbstractRHSCollection* initCollection(std::vector<AbstractMatrixInitializer*> mat_inits,
						std::vector<AbstractCouplingInitializer*> c_inits,
						unsigned int block_size,
						mpi::communicator& comm) {
    return new T(mat_inits, c_inits, block_size, comm);
  }

  double k_x;
  double k_y;
  double k_z;
};


#endif

