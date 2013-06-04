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

  int dummy;
  std::ios_base::openmode om = (world.rank() == 0 && iy == 0) ? std::ios::out : std::ios::app;

  std::ofstream dump("tdtest", om);
  if (world.rank() != 0)
    world.recv(world.rank()-1, 0, dummy);
  for(int i=0; i < 9; i++) {
    dump << test_rhs[0][i] << " ";
  }
  dump.flush();
  if (world.rank()+1 < world.size()) {
    world.send(world.rank()+1, 0, dummy);
  }
  if(world.rank() == 4)
    std::cout << std::endl << std::endl;

  world.barrier();
  crc.doLines(test_rhs);

  for(int il=0; il < 1; il++) {
    if (world.rank() != 0)
      world.recv(world.rank()-1, 0, dummy);
    for(int i=0; i < 9; i++) {
      dump << test_rhs[il][i] << " ";
    }
    dump.flush();
    if (world.rank()+1 < world.size()) {
      world.send(world.rank()+1, 0, dummy);
    }
    // if(world.rank() == 4) {
    //   for(int i=0; i < 9; i++) {
    //   	std::cout << test_rhs[il][i] << " ";
    //   }
    //   std::cout << std::endl << std::endl;
    // }
  }
}
