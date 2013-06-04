#include "RHSCollection.hpp"
#include <fstream>

extern "C" {
  int dgtsv_(int*, int*, double*, double*, double*, double*, int*, int*);
}

int main (int argc, char* argv []) {
  mpi::environment env (argc, argv);
  mpi::communicator world;

  TestMatrixInitializer init;
  std::vector<AbstractMatrixInitializer*> mat_inits(9, & init);
  VacuumCouplingInitializer c_init(& init, 9, world);
  std::vector<AbstractCouplingInitializer*> c_inits(9, & c_init);

  DelegatedRHSCollection crc(mat_inits, c_inits, 9, world);

  double test_rhs_storage[81];
  double* test_rhs[9];

  std::fill_n(test_rhs_storage, 81, 0);

  for(int i = 0; i < 9; i++) {
    test_rhs_storage[i*9 + 0] = 1;
    // test_rhs_storage[i*9 + 8] = 1;
    test_rhs[i] = & test_rhs_storage[i*9];
  }

  world.barrier();
  crc.doLines(test_rhs);

  int dummy;
  std::ios::openmode om = (world.rank() == 0) ? std::ios::out : std::ios::app;

  for(unsigned int il=0; il < 9; il++) {
    std::ofstream dump("tdtest.txt", om);
    if (world.rank() != 0)
      world.recv(world.rank()-1, 0, dummy);
    for(int i=0; i < 9; i++) {
      dump << test_rhs[il][i] << " ";
    }
    dump.flush();
    if (world.rank()+1 < world.size()) {
      world.send(world.rank()+1, 0, dummy);
    }
    if(world.rank() == 4) {
      dump << std::endl << std::endl;
      world.send(0, 0, dummy);
    }
    if(world.rank() == 0)
      world.recv(4, 0, dummy);
    om = std::ios::app;
  }

  if(world.rank() == 4) {
    std::ofstream dump("tdtest.txt", std::ios::app);

    double herp[45];
    double d[45];
    double ld[44];
    double ud[44];

    std::fill_n(herp, 45, 0);
    std::fill_n(d, 45, 1);
    std::fill_n(ld, 44, -1.0/3.0);
    std::fill_n(ud, 44, -1.0/3.0);

    for(unsigned int i=0; i < 45; i += 9) {herp[i] = 1; // herp[i+8] = 1;
    }

    dump << "serial solve: " << std::endl;

    int ninety=45; int one=1; int info;
    dgtsv_(& ninety, & one, ud, d, ld, herp, & ninety, & info);

    for(int i = 0; i < 45; i++)
      dump << herp[i] << " ";
    dump << std::endl;
  }
}
