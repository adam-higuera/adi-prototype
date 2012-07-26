#include "MatrixBlock.hpp"

extern "C" {
  int dgtsv_(int*, int*, double*, double*, double*, double*, int*, int*);
}

int main (int argc, char* argv []) {
  mpi::environment env (argc, argv);
  mpi::communicator world;

  MatrixBlock b(world, 45, ToeplitzMatrixInitializer(1, -1.0/3.0));

  double test_rhs[9];

  std::fill_n(test_rhs, 9, 0);

  test_rhs[0] = 1;
  test_rhs[8] = 1;

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
  if(world.rank() == 4)
    std::cout << std::endl << std::endl;

  world.barrier();
  
  b.solve(test_rhs);

  if (world.rank() != 0)
    world.recv(world.rank()-1, 0, dummy);
  for(int i=0; i < 9; i++) {
    std::cout << test_rhs[i] << " ";
  }
  std::cout.flush();
  if (world.rank()+1 < world.size()) {
    world.send(world.rank()+1, 0, dummy);
  }
  if(world.rank() == 4)
    std::cout << std::endl << std::endl;

  double herp[45];
  double d[45];
  double ld[44];
  double ud[44];

  std::fill_n(herp, 45, 0);
  std::fill_n(d, 45, 1);
  std::fill_n(ld, 44, -1.0/3.0);
  std::fill_n(ud, 44, -1.0/3.0);

  for(unsigned int i=0; i < 45; i += 9) herp[i] = herp[i+8] = 1;

  if(world.rank() == 4) {
    std::cout << "serial solve: " << std::endl;

    int forty_five=45; int one=1; int info;
    dgtsv_(& forty_five, & one, ud, d, ld, herp, & forty_five, & info);

    for(int i = 0; i < 45; i++)
      std::cout << herp[i] << " ";
    std::cout << std::endl;
  }

  double derp[30] = {
    0.99863, 0.987688, 0.965926, 0.93358, 0.891007,
    0.838671, 0.777146, 0.707107, 0.62932, 0.544639,
    0.45399, 0.358368, 0.258819, 0.156434, 0.052336,
    -0.052336, -0.156434, -0.258819, -0.358368, -0.45399,
    -0.544639, -0.62932, -0.707107, -0.777146, -0.838671,
    -0.891007, -0.93358, -0.965926, -0.987688, -0.99863};

  std::fill_n(d, 30, 1 + 2*1.265);
  std::fill_n(ld, 29, -1.265);
  std::fill_n(ud, 29, -1.265);

  double shazbot[29];

  for(unsigned int i=0; i <29; i++) shazbot[i] = -(derp[i+1] - derp[i]) *sqrt(1.265);

  if(world.rank() == 4) {
    std::cout << std::endl << "ARGLE BARGLE: ";
    int twenty_nine=29; int one=1; int info;
    dgtsv_(& twenty_nine, & one, ud, d, ld, shazbot, & twenty_nine, & info);
    for(int i = 0; i < 29; i++)
      std::cout << shazbot[i] << " ";
    std::cout << std::endl;
  }

  std::cout << "Not even doing anything, man!" << std::endl;
}
