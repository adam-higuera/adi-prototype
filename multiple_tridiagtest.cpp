#include "MatrixBlock.hpp"

extern "C" {
  int dgtsv_(int*, int*, double*, double*, double*, double*, int*, int*);
}

int main (int argc, char* argv []) {
  mpi::environment env (argc, argv);
  mpi::communicator world;

  MatrixBlock b(world, 45, ToeplitzMatrixInitializer(1, -1.0/3.0));

  double* test_rhs[9];
  double test_rhs_storage[81];
  for(unsigned int i=0; i < 9; i++)
	test_rhs[i] = & test_rhs_storage[9*i];

  std::fill_n(test_rhs_storage, 27, 0);

  for(int i=0; i < 9; i++) {
	test_rhs[i][0] = 1;
	test_rhs[i][8] = 1;
  }
  world.barrier();
  
  int dummy;
  if (world.rank() != 0)
    world.recv(world.rank()-1, 0, dummy);
  for(int i=0; i < 9; i++) {
    std::cout << test_rhs[0][i] << " ";
  }
  std::cout.flush();
  if (world.rank()+1 < world.size()) {
    world.send(world.rank()+1, 0, dummy);
  }
  if(world.rank() == 4)
    std::cout << std::endl << std::endl;

  world.barrier();

  if(world.rank() == 0)
	std::cout << "Entering solve!" << std::endl;

  b.solveSeveral(test_rhs);

  if(world.rank() == 0)
	std::cout << "Exiting solve!" << std::endl;

  for(int i =0; i <3; i++) {
	world.barrier();

	std::cout << "<<" << world.rank() << ">> ";
	if (world.rank() != 0) 
	  world.recv(world.rank()-1, 0, dummy);
	std::cout << "<<" << world.rank() << ">> ";
	for(int j=0; j < 9; j++) {
	  std::cout << test_rhs[i][j] << " ";
	}
	std::cout.flush();
	if (world.rank()+1 < world.size()) {
	  world.send(world.rank()+1, 0, dummy);
	}
	if(world.rank() == 4)
	  std::cout << std::endl << std::endl;
  }

  if(world.rank() == 4) {
    double herp[45];
    double d[45];
    double ld[45];
    double ud[45];

    std::fill_n(herp, 45, 0);
    std::fill_n(d, 45, 1);
    std::fill_n(ld, 45, -1.0/3.0);
    std::fill_n(ud, 45, -1.0/3.0);

    for(unsigned int i=0; i < 45; i += 9) herp[i] = herp[i+8] = 1;

    std::cout << "serial solve: " << std::endl;

    int ninety=45; int one=1; int info;
    dgtsv_(& ninety, & one, ud, d, ld, herp, & ninety, & info);

    for(int i = 0; i < 45; i++)
      std::cout << herp[i] << " ";
    std::cout << std::endl;
  }
}
