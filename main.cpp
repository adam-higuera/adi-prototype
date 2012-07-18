#include "EMSquare2D.hpp"

int main(int argc, char* argv []) {
  mpi::environment env(argc, argv);
  mpi::communicator world;
  TE10Initializer init(.1/30.0, .1/30.0, .1, .1, 10, world);
  EMSquare2D the_simulation(.1, .1, 5e-10, 30, 2, & init, world);
  the_simulation.simulate();
}
