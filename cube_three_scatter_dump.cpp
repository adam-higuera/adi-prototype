#include "Simulation3D.hpp"
#include <cstdlib>
#include <sstream>
#include <unistd.h>
#include <cmath>

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
  unsigned int n_procs = static_cast<unsigned int>(pow(world.size(), 1.0/3.0));
  unsigned int n_cells = block_size*n_procs;
  ShiftedTEmnlInitializer<1,1,1, ThreeScatterRHSCollection>
    init(.1/n_cells, .1/n_cells, .1/n_cells, .1, .1, .1,
	 block_size);

  double simulation_time=10*sqrt(2)*(.1)/LIGHTSPEED/1.0; // Run for 10 periods
  Simulation3D the_simulation(.1, .1, .1, simulation_time, n_cells, 100,
			      n_procs, n_procs, n_procs, block_size, dump_directory, & init, world);

  the_simulation.simulate(true, 10, 10);
}
