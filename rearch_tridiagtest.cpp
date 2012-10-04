#include "MatrixInitializer.hpp"
#include "CouplingInitializer.hpp"
#include "RHSCollection.hpp"

extern "C" {
  int dgtsv_(int*, int*, double*, double*, double*, double*, int*, int*);
}

int main (int argc, char* argv []) {
  mpi::environment env (argc, argv);
  mpi::communicator world;

  std::vector<AbstractMatrixInitializer*> mat_inits(5, new ToeplitzMatrixInitializer(1, -1.0/3.0));
  std::vector<AbstractCouplingInitializer*> coupling_inits(5, new VacuumCouplingInitializer(mat_inits[0], 5, world));
  CollectiveRHSCollection rhsColl(mat_inits, coupling_inits, 5, world);
  
  double* storage = new double[25];
  double** RHSs = new double*[5];
  for(unsigned int i = 0; i < 5; i++) RHSs[i] = &storage[5*i];
  std::fill_n(storage, 25, 0);
  for(unsigned int i = 0; i < 5; i++) RHSs[i][0] = 1;//RHSs[i][4] = 1;

  rhsColl.doLines(RHSs);

  for(unsigned int ip=0; ip < world.size(); ip++) {
    if (world.rank() == ip)
      for(unsigned int i = 0; i < 5; i++)
  	std::cout << RHSs[0][i] << " ";
    std::cout.flush();
    world.barrier();
  }

  world.barrier();
  world.barrier();
  world.barrier();
  world.barrier();

  if(world.rank() == 4)
    std::cout << std::endl;


  if(world.rank() == 4) {
    double herp[25];
    double d[25];
    double ld[25];
    double ud[25];

    std::fill_n(herp, 25, 0);
    std::fill_n(d, 25, 1);
    std::fill_n(ld, 24, -1.0/3.0);
    std::fill_n(ud, 24, -1.0/3.0);

    for(unsigned int i=0; i < 25; i += 5) herp[i] = 1;// herp[i+4] = 1;

    std::cout << "serial solve: " << std::endl;

    int ninety=25; int one=1; int info;
    dgtsv_(& ninety, & one, ud, d, ld, herp, & ninety, & info);

    for(int i = 0; i < 25; i++)
      std::cout << herp[i] << " ";
    std::cout << std::endl;
  }
}
