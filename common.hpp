#ifndef MY_COMMON_H
#define MY_COMMON_H

#include <boost/mpi.hpp>
#define LIGHTSPEED 2.977e8

namespace mpi = boost::mpi;

extern "C" {
  int dgttrf_(int*, double*, double*, double*, double*, int*, int*);
  int dgttrs_(char*, int*, int*, double*, double*, double*, double*, int*, double*, int*, int*);
}

enum WhichDiagonal {MAIN_DIAG, UPPER_DIAG, LOWER_DIAG};
enum BoundaryLocation {NO_BDY, UPPER_BDY, LOWER_BDY, DOUBLE_BDY};
enum bdyDir {BDY_X, BDY_Y};

#endif
