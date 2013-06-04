#ifndef MATRIX_INITIALIZER_H
#define MATRIX_INITIALIZER_H

#include "common.hpp"

class AbstractMatrixInitializer {
public:
  AbstractMatrixInitializer() {}
  ~AbstractMatrixInitializer() {}

  virtual double operator()(unsigned int i, WhichDiagonal w)=0;
};

class ToeplitzMatrixInitializer : public AbstractMatrixInitializer {
public:
  ToeplitzMatrixInitializer(double diag, double off_diag) : diag (diag), offDiag(off_diag) {}
  double operator ()(unsigned int i, WhichDiagonal w) {
	switch (w) {
	case MAIN_DIAG:
	  return this->diag;
	case LOWER_DIAG:
	case UPPER_DIAG:
	  return this->offDiag;
	}
	return -1.0;
  }
protected:
  double diag, offDiag;
};

class VacuumMatrixInitializer : public ToeplitzMatrixInitializer {
public:
  VacuumMatrixInitializer(double dx, double dt, unsigned int block_size, BoundaryLocation bdy_loc) :
    ToeplitzMatrixInitializer(1 + 2.0*LIGHTSPEED*LIGHTSPEED*dt*dt/((2*dx)*(2*dx)),
			      -LIGHTSPEED*LIGHTSPEED*dt*dt/((2*dx)*(2*dx))),
    blockSize(block_size),
    bdyLoc(bdy_loc) {}

  double operator ()(unsigned int i, WhichDiagonal w) {
    if (this->checkBdy(i, w))
      return this->diag + this->offDiag;
    else
      return ToeplitzMatrixInitializer::operator()(i, w);
    return -1.0;
  }

private:
  bool checkBdy(unsigned int i, WhichDiagonal w) {
    return (w == MAIN_DIAG &&
	    (i == 0 && (bdyLoc == LOWER_BDY || bdyLoc == DOUBLE_BDY)) ||
	    (i == blockSize - 1 && (bdyLoc == UPPER_BDY || bdyLoc == DOUBLE_BDY)));
  }
  BoundaryLocation bdyLoc;
  unsigned int blockSize;
};

class TestMatrixInitializer : public AbstractMatrixInitializer {
public:
  double operator ()(unsigned int i, WhichDiagonal w) {
    switch(w) {
    case MAIN_DIAG:
      return 1.0;
    case LOWER_DIAG:
    case UPPER_DIAG:
      return -1.0/3.0;
    }
    return -1.0;
  }
};

#endif
