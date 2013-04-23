all: one-line-yee one-line-local-only one-line-comm-only one-line-reduced-only one-line-full one-line-chunked-local-only one-line-chunked-comm-only one-line-chunked-reduced-only one-line-chunked-full \
square-local-only square-comm-only square-reduced-only square-full square-chunked-local-only square-chunked-comm-only square-chunked-reduced-only square-chunked-full square-full-dump square-reduced-only-dump square-comm-only-dump square-local-only-dump

one-line-yee: one_line.cpp
	mpicxx -g -DYEE -o $@ -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include one_line.cpp Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp -lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

one-line-local-only: one_line.cpp
	mpicxx -g -DNO_NEAREST_NEIGHBOR -DNO_COLLECTIVES -DNO_REDUCED_SOLVE -o $@ -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include one_line.cpp Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp -lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

one-line-comm-only: one_line.cpp
	mpicxx -g -DCOMMUNICATION_ONLY -DNO_LOCAL_SOLVES -DNO_REDUCED_SOLVE -o $@ -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include one_line.cpp Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp -lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

one-line-reduced-only: one_line.cpp
	mpicxx -g -DNO_NEAREST_NEIGHBOR -DNO_COLLECTIVES -DNO_EXPLICIT_SOLVE -DNO_LOCAL_SOLVES -o $@ -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include one_line.cpp Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp -lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

one-line-full: one_line.cpp
	mpicxx -g -o $@ -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include one_line.cpp Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp -lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

one-line-chunked-local-only: one_line_chunked.cpp
	mpicxx -g -DNO_NEAREST_NEIGHBOR -DNO_COLLECTIVES -DNO_REDUCED_SOLVE -o $@ -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include one_line_chunked.cpp Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp -lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

one-line-chunked-comm-only: one_line_chunked.cpp
	mpicxx -g -DCOMMUNICATION_ONLY -DNO_LOCAL_SOLVES -DNO_REDUCED_SOLVE -o $@ -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include one_line_chunked.cpp Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp -lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

one-line-chunked-reduced-only: one_line_chunked.cpp
	mpicxx -g -DNO_NEAREST_NEIGHBOR -DNO_COLLECTIVES -DNO_EXPLICIT_SOLVE -DNO_LOCAL_SOLVES -o $@ -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include one_line_chunked.cpp Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp -lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

one-line-chunked-full: one_line_chunked.cpp
	mpicxx -g -o $@ -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include one_line_chunked.cpp Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp -lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

square-local-only: square_collective.cpp
	mpicxx -g -DNO_NEAREST_NEIGHBOR -DNO_COLLECTIVES -DNO_REDUCED_SOLVE -o $@ -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include square_collective.cpp Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp -lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

square-comm-only: square_collective.cpp
	mpicxx -g -DCOMMUNICATION_ONLY -DNO_LOCAL_SOLVES -DNO_REDUCED_SOLVE -o $@ -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include square_collective.cpp Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp -lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

square-reduced-only: square_collective.cpp
	mpicxx -g -DNO_NEAREST_NEIGHBOR -DNO_COLLECTIVES -DNO_EXPLICIT_SOLVE -DNO_LOCAL_SOLVES -o $@ -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include square_collective.cpp Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp -lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

square-full: square_collective.cpp
	mpicxx -g -o $@ -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include square_collective.cpp Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp -lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

square-chunked-local-only: square_chunked.cpp
	mpicxx -g -DNO_NEAREST_NEIGHBOR -DNO_COLLECTIVES -DNO_REDUCED_SOLVE -o $@ -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include square_chunked.cpp Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp -lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

square-chunked-comm-only: square_chunked.cpp
	mpicxx -g -DCOMMUNICATION_ONLY -DNO_LOCAL_SOLVES -DNO_REDUCED_SOLVE -o $@ -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include square_chunked.cpp Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp -lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

square-chunked-reduced-only: square_chunked.cpp
	mpicxx -g -DNO_NEAREST_NEIGHBOR -DNO_COLLECTIVES -DNO_EXPLICIT_SOLVE -DNO_LOCAL_SOLVES -o $@ -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include square_chunked.cpp Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp -lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

square-chunked-full: square_chunked.cpp
	mpicxx -g -o $@ -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include square_chunked.cpp Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp -lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

square-local-only-dump: square_dump.cpp
	mpicxx -g -DNO_NEAREST_NEIGHBOR -DNO_COLLECTIVES -DNO_REDUCED_SOLVE -o $@ -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include square_dump.cpp Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp -lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

square-comm-only-dump: square_dump.cpp
	mpicxx -g -DCOMMUNICATION_ONLY -DNO_LOCAL_SOLVES -DNO_REDUCED_SOLVE -o $@ -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include square_dump.cpp Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp -lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

square-reduced-only-dump: square_dump.cpp
	mpicxx -g -DNO_NEAREST_NEIGHBOR -DNO_COLLECTIVES -DNO_EXPLICIT_SOLVE -DNO_LOCAL_SOLVES -o $@ -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include square_dump.cpp Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp -lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

square-full-dump: square_dump.cpp
	mpicxx -g -o $@ -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include square_dump.cpp Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp -lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread