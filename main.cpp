#include "EMSquare2D.hpp"

int main(int argc, char* argv []) {
  mpi::environment env(argc, argv);
  mpi::communicator world;
  unsigned int simulation_size = 120.0*7*3*2; // LCM of 2,3,4,5,6,7,8,9
  unsigned int block_size = static_cast<unsigned int>(simulation_size/sqrt(world.size()));
  TEmnInitializer<1,1> init(.1/(simulation_size), .1/(simulation_size), .1, .1, block_size, world);
  EMSquare2D the_simulation(.1, .1, 5e-8, 120*7, 1000, & init, world);
  the_simulation.simulate(false);
}
