#include "MatrixBlock.hpp"

int main (int argc, char* argv []) {
  mpi::environment env (argc, argv);
  mpi::communicator world;

  MatrixBlock b(world, 45, ToeplitzMatrixInitializer(1, -1.0/3.0));

  double test_rhs[9];

  std::fill_n(test_rhs, 9, 0);

  test_rhs[5] = 5 + world.rank()*9;

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

  for(unsigned int i=5; i < 45; i += 9) herp[i] = i;

  if(world.rank() == 4) {
	std::cout << "serial solve: " << std::endl;

	int forty_five=45; int one=1; int info;
	dgtsv_(& forty_five, & one, ud, d, ld, herp, & forty_five, & info);

	for(int i = 0; i < 45; i++)
	  std::cout << herp[i] << " ";
	std::cout << std::endl;
  }
}
