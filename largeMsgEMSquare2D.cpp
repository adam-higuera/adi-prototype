#include "largeMsgEMSquare2D.hpp"

largeMsgEMSquare2D::largeMsgEMSquare2D(
		       double L_x, double L_y, double T,
		       unsigned int n_cells, unsigned int n_steps,
		       EMSquare2DInitializer* init, mpi::communicator & world
		       )
: EMSquare2D(L_x, L_y, T, n_cells, n_steps, init, world),
  implicit_rhs_holder(new double*[this->block_size])
{
  implicit_rhs_holder[0] = new double[this->block_size*this->block_size];
  for(unsigned int il=0; il < this->block_size - 1; il++)
    implicit_rhs_holder[il+1] = implicit_rhs_holder[il] + this->block_size;
}

void largeMsgEMSquare2D::implicitUpdateM() {
  for(unsigned int iy=0; iy < this->block_size; iy++) {
    for(unsigned int ix=0; ix < this->block_size; ix++) {
      implicit_rhs_holder[iy][ix] = B_z[ix][iy] - LIGHTSPEED*dt/(2*dy) * (E_y[ix+1][iy] - E_y[ix][iy]);
    }
  }

  Implicit_dx.solveSeveral(implicit_rhs_holder);

  for(unsigned int iy=0; iy < this->block_size; iy++) {
    for(unsigned int ix=0; ix < this->block_size; ix++) {
      B_z[ix][iy] = implicit_rhs_holder[iy][ix];
    }
  }

  this->implicitMSubstituteB();
}

void largeMsgEMSquare2D::implicitUpdateP() {
  for(unsigned int ix=0; ix < this->block_size; ix++) {
    for(unsigned int iy=0; iy < this->block_size; iy++) {
      implicit_rhs_holder[ix][iy] = B_z[ix][iy] + LIGHTSPEED*dt/(2*dy) * (E_x[ix][iy+1] - E_x[ix][iy]);
    }
  }

  Implicit_dy.solveSeveral(implicit_rhs_holder);

  std::copy(implicit_rhs_holder[0],  implicit_rhs_holder[0] + block_size*block_size, B_z[0]);

  this->implicitPSubstituteB();
}
