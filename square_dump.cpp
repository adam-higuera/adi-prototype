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

  unsigned int block_size = domain_size ? std::atoi(domain_size) : 50;
  std::string dump_directory = dump_dir ? std::string(dump_dir) : std::string("");
  if (world.rank() == 0)
    std::cout << "domain_size: " << block_size << std::endl;
  unsigned int n_procs = static_cast<unsigned int>(sqrt(world.size()));
  unsigned int n_cells = block_size*n_procs;
  TEmnInitializer<1,1, CollectiveRHSCollection>
    init(.1/n_cells, .1/n_cells, .1, .1, n_procs, n_procs, block_size, world);
  Simulation the_simulation(.1, .1, 5e-8, n_cells, 10000, n_procs, n_procs, block_size, dump_directory, & init, world);
  the_simulation.simulate(true, 1, 10);
}
