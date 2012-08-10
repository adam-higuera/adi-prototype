#include "MatrixBlock.hpp"

extern "C" {
  int dgtsv_(int*, int*, double*, double*, double*, double*, int*, int*);
}

int main (int argc, char* argv []) {
  mpi::environment env (argc, argv);
  mpi::communicator world;

  MatrixBlock b(world, 225, ToeplitzMatrixInitializer(1, -1.0/3.0));

  double test_rhs[9];

  std::fill_n(test_rhs, 9, 0);

  test_rhs[0] = 1;
  test_rhs[8] = 1;

  world.barrier();
  
  int dummy;
  if (world.rank() != 0)
    world.recv(world.rank()-1, 0, dummy);
  for(int i=0; i < 9; i++) {
    std::cout << test_rhs[i] << " ";
  }
  std::cout.flush();
  if (world.rank()+1 < world.size()) {
    world.send(world.rank()+1, 0, dummy);
  }
  if(world.rank() == 24)
    std::cout << std::endl << std::endl;

  world.barrier();

  b.solve(test_rhs);

  world.barrier();

  std::cout << "<<" << world.rank() << ">> ";
  if (world.rank() != 0) 
    world.recv(world.rank()-1, 0, dummy);
  std::cout << "<<" << world.rank() << ">> ";
  for(int i=0; i < 9; i++) {
    std::cout << test_rhs[i] << " ";
  }
  std::cout.flush();
  if (world.rank()+1 < world.size()) {
    world.send(world.rank()+1, 0, dummy);
  }
  if(world.rank() == 24)
    std::cout << std::endl << std::endl;

  if(world.rank() == 24) {
    double herp[225];
    double d[225];
    double ld[224];
    double ud[224];

    std::fill_n(herp, 225, 0);
    std::fill_n(d, 225, 1);
    std::fill_n(ld, 224, -1.0/3.0);
    std::fill_n(ud, 224, -1.0/3.0);

    for(unsigned int i=0; i < 225; i += 9) herp[i] = herp[i+8] = 1;

    std::cout << "serial solve: " << std::endl;

    int ninety=225; int one=1; int info;
    dgtsv_(& ninety, & one, ud, d, ld, herp, & ninety, & info);

    for(int i = 0; i < 225; i++)
      std::cout << herp[i] << " ";
    std::cout << std::endl;
  }
}
