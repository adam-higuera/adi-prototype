include "RHSCollection.hpp"
#include <fstream>
#include <iomanip>

extern "C" {
  int dgtsv_(int*, int*, double*, double*, double*, double*, int*, int*);
}

int main (int argc, char* argv []) {
  mpi::environment env (argc, argv);
  mpi::communicator world;

  TestMatrixInitializer init;
  std::vector<AbstractMatrixInitializer*> mat_inits(81, & init);
  VacuumCouplingInitializer c_init(& init, 81, world);
  std::vector<AbstractCouplingInitializer*> c_inits(81, & c_init);

  std::cout << "Entering constructor" << world.rank() << std::endl;
  CollectiveRHSCollection crc(mat_inits, c_inits, 9, world);
  std::cout << "Exiting constructor" << world.rank() << std::endl;

  double test_rhs_storage[81*9];
  double* test_rhs[81];

  std::fill_n(test_rhs_storage, 81*9, 0);

  for(int i = 0; i < 81; i++) {
    test_rhs_storage[i*9 + 0] = 1;
    // test_rhs_storage[i*9 + 8] = 1;
    test_rhs[i] = & test_rhs_storage[i*9];
  }

  int dummy;
  // std::ios_base::openmode om = (world.rank() == 0 && iy == 0) ? std::ios::out : std::ios::app;

  // std::ofstream dump("tdtest", om);
  // if (world.rank() != 0)
  //   world.recv(world.rank()-1, 0, dummy);
  // for(int i=0; i < 9; i++) {
  //   dump << test_rhs[0][i] << " ";
  // }
  // dump.flush();
  // if (world.rank()+1 < world.size()) {
  //   world.send(world.rank()+1, 0, dummy);
  // }
  // if(world.rank() == 4)
  //   std::cout << std::endl << std::endl;

  world.barrier();
  std::cout << "Entering do_lines" << world.rank() << std::endl;
  crc.doLines(test_rhs);
  std::cout << "Exiting do_lines" << world.rank() << std::endl;

  for(unsigned int il = 0; il < 81; il++) {
    if(world.rank() != 0)
      world.recv(world.rank() - 1, 0, dummy);
    else if (il != 0)
      world.recv(world.size() - 1, 0, dummy);

    std::ios_base::openmode om = (world.rank() == 0 && il == 0) ? std::ios::out : std::ios::app;

    std::ofstream dump("tdtest.txt", om);
    for(unsigned int ix = 0; ix < 9; ix++) {
      dump << test_rhs[il][ix];
      if (ix != 9-1 || world.rank() != world.size() - 1)
	dump << " ";
    }
    if(world.rank() == world.size() - 1)
      dump << "\n";
    dump.flush();
    if(world.rank() != world.size() - 1)
      world.send(world.rank()+1, 0, dummy);
    else if (il != 81 - 1)
      world.send(0, 0, dummy);
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
