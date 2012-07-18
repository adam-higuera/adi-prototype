#include "MatrixBlock.hpp"
#include <boost/mpi.hpp>
#define C 2.977e8

enum bdy_dir {BDY_X, BDY_Y};

class EMSquare2DInitializer {
public:
  EMSquare2DInitializer(double dx, double dy, double L_x, double L_y,
						unsigned int block_size, const mpi::communicator & comm);
  virtual double E_x(unsigned int i, unsigned int j) = 0;
  virtual double E_y(unsigned int i, unsigned int j) = 0;
  virtual double B_z(unsigned int i, unsigned int j) = 0;
protected:
  double x(unsigned int i);
  double y(unsigned int j);
  
  double dx;
  double dy;
  double x_offset;
  double y_offset;
  double L_x;
  double L_y;
};

class EMSquare2D {
public:
  EMSquare2D(double L_x, double L_y, double T,
			 unsigned int n_cells, unsigned int n_steps,
			 EMSquare2DInitializer* init, mpi::communicator & world);

  void simulate();

private:
  mpi::communicator & world, x_line, y_line;

  unsigned int block_size, n_steps;

  double** E_x;
  double** E_y;
  double** B_z;

  double* rhs_holder;

  double dx;
  double dy;
  double dt;

  double* boundary_E_in;
  double* boundary_B_in;
  double* boundary_E_out;
  double* boundary_B_out;
  
  MatrixBlock Implicit_dy;
  MatrixBlock Implicit_dx;

  double get_bdy_B_z(mpi::communicator comm);
  void exchange_bdy_values(mpi::communicator comm, bdy_dir dir);

  void implicitUpdateM();
  void explicitUpdateP();
  void implicitUpdateP();
  void explicitUpdateM();

  void TimeStep();

  void dumpFields(std::string filename);
};

class TE10Initializer : public EMSquare2DInitializer {
public:
  TE10Initializer(double dx, double dy, double L_x, double L_y,
				  unsigned int block_size, const mpi::communicator & comm);
  virtual double E_x(unsigned int i, unsigned int j);
  virtual double E_y(unsigned int i, unsigned int j);
  virtual double B_z(unsigned int i, unsigned int j);
};
