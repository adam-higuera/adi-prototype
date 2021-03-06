#include "Simulation.hpp"
#include <cstdlib>
#include <sstream>
#include <unistd.h>

int main(int argc, char* argv []) {  
  mpi::environment env(argc, argv);
  mpi::communicator world;

  // std::ostringstream startstrm2;
  // char hostname[520];
  // gethostname(hostname, 512);
  // startstrm2 << "Rank " << world.rank() << " is running on host: " <<
  //   hostname << std::endl;
  // std::cout << startstrm2.str() << std::endl;

  char *domain_size = std::getenv("DOMAIN_SIZE");
  char *dump_dir = std::getenv("RESULTS_DIR");

  unsigned int block_size = domain_size ? std::atoi(domain_size) : 50;
  std::string dump_directory = dump_dir ? std::string(dump_dir) : std::string("");

  if (world.rank() == 0)
    std::cout << "domain_size: " << block_size << std::endl;
  unsigned int n_cells = block_size*world.size();
  TEmnInitializer<1,1, ChunkedRHSCollection> init(.1/n_cells, .1/block_size, .1, .1, world.size(), 1, block_size, world);
  /*
		       double L_x, double L_y, double T,
		       unsigned int n_cells, unsigned int n_steps,
		       unsigned int procs_x, unsigned int procs_y,
		       SimulationInitializer* init, mpi::communicator & world
  */
  double simulation_time=10*sqrt(2)*(.1)/LIGHTSPEED; // Run for 100 periods
  Simulation the_simulation(.1, .1, simulation_time, n_cells, 10000,
			    world.size(), 1, block_size, dump_directory, & init, world);
  the_simulation.simulate(false, 100, 10);
}
