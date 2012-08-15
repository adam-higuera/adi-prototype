#include "largeMsgEMSquare2D.hpp"
#include <cstdlib>
#include <sstream>
#include <unistd.h>
#include <boost/program_options.hpp>

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

  unsigned int block_size = domain_size ? std::atoi(domain_size) : 1000;
  if (world.rank() == 0)
    std::cout << "domain_size: " << block_size << std::endl;
  unsigned int n_cells = block_size*static_cast<unsigned int>(sqrt(world.size()));
  TEmnInitializer<1,1> init(.1/n_cells, .1/n_cells, .1, .1, block_size, world);
  largeMsgEMSquare2D the_simulation(.1, .1, 5e-8, n_cells, 1000, & init, world);
  the_simulation.simulate(false);
}
