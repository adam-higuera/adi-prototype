#include "MatrixBlock.hpp"
#include <boost/mpi.hpp>
#include "EMSquare2D.hpp"

class largeMsgEMSquare2D : public EMSquare2D {
public:
  largeMsgEMSquare2D(double L_x, double L_y, double T,
	     unsigned int n_cells, unsigned int n_steps,
	     EMSquare2DInitializer* init, mpi::communicator & world);

protected:
  virtual void implicitUpdateM();
  virtual void implicitUpdateP();

private:
  double** implicit_rhs_holder;
};
