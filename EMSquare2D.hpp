#include "MatrixBlock.hpp"
#include <boost/mpi.hpp>
#define LIGHTSPEED 2.977e8

enum bdy_dir {BDY_X, BDY_Y};

class EMSquare2DInitializer {
public:
  EMSquare2DInitializer(double dx, double dy, double L_x, double L_y,
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

class EMSquare2D {
public:
  EMSquare2D(double L_x, double L_y, double T,
	     unsigned int n_cells, unsigned int n_steps,
	     EMSquare2DInitializer* init, mpi::communicator & world);

  void simulate(bool dump=true, unsigned int dump_periodicity=9, unsigned int total_dumps=100);

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

  double* B_left_bdy;
  double* B_right_bdy;
  double* B_bot_bdy;
  double* B_top_bdy;

  double* bdy_out_buffer_topright;
  double* bdy_out_buffer_botleft;
  
  MatrixBlock Implicit_dy;
  MatrixBlock Implicit_dx;

  double get_bdy_B_z(mpi::communicator comm);
  void exchange_bdy_values(mpi::communicator comm, bdy_dir dir);

  void implicitUpdateM();
  void explicitUpdateP();
  void implicitUpdateP();
  void explicitUpdateM();

  void TimeStep();
  void printField(std::string msg);

  void dumpFields(std::string filename);
};

template<int m, int n>
class TEmnInitializer : public EMSquare2DInitializer {
public:
  TEmnInitializer(double dx, double dy, double L_x, double L_y,
		  unsigned int block_size, const mpi::communicator & comm)
    : EMSquare2DInitializer(dx, dy, L_x, L_y, block_size, comm) {}
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
