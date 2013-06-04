LIBS=-lboost_mpi -lboost_program_options -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread
LIB_DIR=-L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib
INC_DIR=-I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include
SRC_FILES=Simulation.cpp RHSCollection.cpp LocalReducedRHS.cpp RemoteReducedRHS.cpp TDCoupling.cpp LocalSolver.cpp
LOCAL_ONLY_FLAGS=-DNO_NEAREST_NEIGHBOR -DNO_COLLECTIVES -DNO_REDUCED_SOLVE
COMM_ONLY_FLAGS=-DCOMMUNICATION_ONLY -DNO_LOCAL_SOLVES -DNO_REDUCED_SOLVE
REDUCED_ONLY_FLAGS=-DCOMMUNICATION_ONLY -DNO_NEAREST_NEIGHBOR -DNO_COLLECTIVES -DNO_EXPLICIT_SOLVE -DNO_LOCAL_SOLVES
TOTAL_REDUCED_ONLY_FLAGS=-DTOTAL_REDUCED_ONLY -DNO_EXPLICIT_SOLVE
EVERYTHING_ELSE_FLAGS=-DNO_TOTAL_REDUCED
HEADERS=Simulation.hpp RHSCollection.hpp

all: one-line-yee one-line-local-only one-line-comm-only one-line-reduced-only one-line-full \
one-line-chunked-local-only one-line-chunked-comm-only one-line-chunked-reduced-only one-line-chunked-full \
square-local-only square-comm-only square-reduced-only square-full square-chunked-local-only \
square-chunked-comm-only square-chunked-reduced-only square-chunked-full square-full-dump \
square-reduced-only-dump square-comm-only-dump square-local-only-dump \
one-line-local-only-barrier one-line-comm-only-barrier one-line-reduced-only-barrier one-line-full-barrier \
one-line-chunked-local-only-barrier one-line-chunked-comm-only-barrier one-line-chunked-reduced-only-barrier one-line-chunked-full-barrier \
square-local-only-barrier square-comm-only-barrier square-reduced-only-barrier square-full-barrier square-chunked-local-only-barrier \
square-chunked-comm-only-barrier square-chunked-reduced-only-barrier square-chunked-full-barrier

one-line-yee: one_line.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -DYEE -o $@ ${LIB_DIR} ${INC_DIR} one_line.cpp ${SRC_FILES} ${LIBS}

one-line-local-only: one_line.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g ${LOCAL_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} one_line.cpp ${SRC_FILES} ${LIBS}

one-line-local-only-barrier: one_line.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -DUSE_BARRIERS ${LOCAL_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} one_line.cpp ${SRC_FILES} ${LIBS}

one-line-comm-only: one_line.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g ${COMM_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} one_line.cpp ${SRC_FILES} ${LIBS}

one-line-comm-only-barrier: one_line.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -DUSE_BARRIERS ${COMM_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} one_line.cpp ${SRC_FILES} ${LIBS}

one-line-reduced-only: one_line.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g ${REDUCED_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} one_line.cpp ${SRC_FILES} ${LIBS}

one-line-reduced-only-barrier: one_line.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -DUSE_BARRIERS ${REDUCED_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} one_line.cpp ${SRC_FILES} ${LIBS}

one-line-full: one_line.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -o $@ ${LIB_DIR} ${INC_DIR} one_line.cpp ${SRC_FILES} ${LIBS}

one-line-full-barrier: one_line.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -DUSE_BARRIERS -o $@ ${LIB_DIR} ${INC_DIR} one_line.cpp ${SRC_FILES} ${LIBS}

one-line-chunked-local-only: one_line_chunked.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g ${LOCAL_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} one_line_chunked.cpp ${SRC_FILES} ${LIBS}

one-line-chunked-local-only-barrier: one_line_chunked.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -DUSE_BARRIERS ${LOCAL_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} one_line_chunked.cpp ${SRC_FILES} ${LIBS}

one-line-chunked-comm-only: one_line_chunked.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g ${COMM_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} one_line_chunked.cpp ${SRC_FILES} ${LIBS}

one-line-chunked-comm-only-barrier: one_line_chunked.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -DUSE_BARRIERS ${COMM_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} one_line_chunked.cpp ${SRC_FILES} ${LIBS}

one-line-chunked-reduced-only: one_line_chunked.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g ${REDUCED_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} one_line_chunked.cpp ${SRC_FILES} ${LIBS}

one-line-chunked-reduced-only-barrier: one_line_chunked.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -DUSE_BARRIERS ${REDUCED_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} one_line_chunked.cpp ${SRC_FILES} ${LIBS}

one-line-chunked-full: one_line_chunked.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -o $@ ${LIB_DIR} ${INC_DIR} one_line_chunked.cpp ${SRC_FILES} ${LIBS}

one-line-chunked-full-barrier: one_line_chunked.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -DUSE_BARRIERS -o $@ ${LIB_DIR} ${INC_DIR} one_line_chunked.cpp ${SRC_FILES} ${LIBS}

square-local-only: square_collective.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g ${LOCAL_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} square_collective.cpp ${SRC_FILES} ${LIBS}

square-local-only-barrier: square_collective.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -DUSE_BARRIERS ${LOCAL_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} square_collective.cpp ${SRC_FILES} ${LIBS}

square-comm-only: square_collective.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g ${COMM_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} square_collective.cpp ${SRC_FILES} ${LIBS}

square-comm-only-barrier: square_collective.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -DUSE_BARRIERS ${COMM_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} square_collective.cpp ${SRC_FILES} ${LIBS}

square-reduced-only: square_collective.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g ${REDUCED_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} square_collective.cpp ${SRC_FILES} ${LIBS}

square-reduced-only-barrier: square_collective.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -DUSE_BARRIERS ${REDUCED_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} square_collective.cpp ${SRC_FILES} ${LIBS}

square-full: square_collective.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -o $@ ${LIB_DIR} ${INC_DIR} square_collective.cpp ${SRC_FILES} ${LIBS}

square-full-barrier: square_collective.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -DUSE_BARRIERS -o $@ ${LIB_DIR} ${INC_DIR} square_collective.cpp ${SRC_FILES} ${LIBS}

square-chunked-local-only: square_chunked.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g ${LOCAL_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} square_chunked.cpp ${SRC_FILES} ${LIBS}

square-chunked-local-only-barrier: square_chunked.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -DUSE_BARRIERS ${LOCAL_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} square_chunked.cpp ${SRC_FILES} ${LIBS}

square-chunked-comm-only: square_chunked.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g ${COMM_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} square_chunked.cpp ${SRC_FILES} ${LIBS}

square-chunked-comm-only-barrier: square_chunked.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -DUSE_BARRIERS ${COMM_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} square_chunked.cpp ${SRC_FILES} ${LIBS}

square-chunked-reduced-only: square_chunked.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g ${REDUCED_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} square_chunked.cpp ${SRC_FILES} ${LIBS}

square-chunked-reduced-only-barrier: square_chunked.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -DUSE_BARRIERS ${REDUCED_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} square_chunked.cpp ${SRC_FILES} ${LIBS}

square-chunked-full: square_chunked.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -o $@ ${LIB_DIR} ${INC_DIR} square_chunked.cpp ${SRC_FILES} ${LIBS}

square-chunked-full-barrier: square_chunked.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -DUSE_BARRIERS -o $@ ${LIB_DIR} ${INC_DIR} square_chunked.cpp ${SRC_FILES} ${LIBS}

square-total-reduced-only: square_collective.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -o $@ ${TOTAL_REDUCED_ONLY_FLAGS} ${LIB_DIR} ${INC_DIR} square_collective.cpp ${SRC_FILES} ${LIBS}

square-everything-else: square_collective.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -o $@ ${EVERYTHING_ELSE_FLAGS} ${LIB_DIR} ${INC_DIR} square_collective.cpp ${SRC_FILES} ${LIBS}

square-chunked-total-reduced-only: square_chunked.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -o $@ ${TOTAL_REDUCED_ONLY_FLAGS} ${LIB_DIR} ${INC_DIR} square_chunked.cpp ${SRC_FILES} ${LIBS}

square-chunked-everything-else: square_chunked.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -o $@ ${EVERYTHING_ELSE_FLAGS} ${LIB_DIR} ${INC_DIR} square_chunked.cpp ${SRC_FILES} ${LIBS}

one-line-total-reduced-only: one_line.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -o $@ ${TOTAL_REDUCED_ONLY_FLAGS} ${LIB_DIR} ${INC_DIR} one_line.cpp ${SRC_FILES} ${LIBS}

one-line-everything-else: one_line.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -o $@ ${EVERYTHING_ELSE_FLAGS} ${LIB_DIR} ${INC_DIR} one_line.cpp ${SRC_FILES} ${LIBS}

one-line-chunked-total-reduced-only: one_line_chunked.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -o $@ ${TOTAL_REDUCED_ONLY_FLAGS} ${LIB_DIR} ${INC_DIR} one_line_chunked.cpp ${SRC_FILES} ${LIBS}

one-line-chunked-everything-else: one_line_chunked.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -o $@ ${EVERYTHING_ELSE_FLAGS} ${LIB_DIR} ${INC_DIR} one_line_chunked.cpp ${SRC_FILES} ${LIBS}

square-local-only-dump: square_dump.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g ${LOCAL_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} square_dump.cpp ${SRC_FILES} ${LIBS}

square-comm-only-dump: square_dump.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g ${COMM_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} square_dump.cpp ${SRC_FILES} ${LIBS}

square-reduced-only-dump: square_dump.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g ${REDUCED_ONLY_FLAGS} -o $@ ${LIB_DIR} ${INC_DIR} square_dump.cpp ${SRC_FILES} ${LIBS}

square-full-dump: square_dump.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -o $@ ${LIB_DIR} ${INC_DIR} square_dump.cpp ${SRC_FILES} ${LIBS}

accuracy-check: accuracy_check.cpp ${HEADERS} ${SRC_FILES}
	mpicxx -g -o $@ ${LIB_DIR} ${INC_DIR} accuracy_check.cpp ${SRC_FILES} ${LIBS}