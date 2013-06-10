#include "Simulation.hpp"
#include <cstdlib>
#include <sstream>
#include <string>
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
  char *n_steps_env = std::getenv("NUM_STEPS");

  unsigned int block_size = domain_size ? std::atoi(domain_size) : 50;
  std::string dump_directory = dump_dir ? std::string(dump_dir) : std::string("");
  unsigned int n_steps = n_steps_env ? std::atoi(n_steps_env) : 10000;
  if (world.rank() == 0)
    std::cout << "domain_size: " << block_size << std::endl;
  unsigned int n_procs = static_cast<unsigned int>(sqrt(world.size()));
  unsigned int n_cells = block_size*n_procs;
  TEmnInitializer<1,1, DelegatedRHSCollection>
    init(.1/n_cells, .1/n_cells, .1, .1, n_procs, n_procs, block_size, world);

  double simulation_time=10*sqrt(2)*(.1)/LIGHTSPEED; // Run for 10 periods
  Simulation the_simulation(.1, .1, simulation_time, n_cells, n_steps,
			    n_procs, n_procs, block_size, dump_directory, & init, world);
  the_simulation.simulate(true, n_steps, 10);
}
