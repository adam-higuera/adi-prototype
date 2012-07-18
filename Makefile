MPI_LIB=mpi
MPICXX_LIB=mpi_cxx
BOOST_MPI_LIB=boost_mpi

HOME=/Users/adamvh
BOOST_ROOT_DIR=${HOME}/Library/boost_1_50_0
BOOST_LIB_DIR=${BOOST_ROOT_DIR}/stage/lib

MPICXX=mpicxx

CXXFLAGS=-I${BOOST_ROOT_DIR}

LDFLAGS=-L${BOOST_LIB_DIR} -l${MPI_LIB} -l${MPICXX_LIB}  -l${BOOST_MPI_LIB} -framework Accelerate

tridiag: tridiag.cpp
	${MPICXX} ${CXXFLAGS} -o $@ tridiag.cpp

tridiagtest: tridiagtest.cpp MatrixBlock.o
	${MPICXX} $^ ${CXXFLAGS} ${LDFLAGS} -o $@

emsquare2d: EMSquare2D.o MatrixBlock.o main.cpp
	${MPICXX} $^ ${CXXFLAGS} ${LDFLAGS} -o $@

MatrixBlock.o: MatrixBlock.cpp MatrixBlock.hpp
	${MPICXX} $< -c ${CXXFLAGS}  -o $@

EMSquare2D.o: EMSquare2D.cpp MatrixBlock.o EMSquare2D.hpp
	${MPICXX} $< -c ${CXXFLAGS}  -o $@

# all: bin/dta2json bin/import_dta_private.mex${MEXSUFFIX}

# bin/dta2json: ${COMMON_OBJS} ${YAJL_OBJS} json/main.o json/dta_to_json.o
# 	${CXX} $^ ${CXXFLAGS} -o $@

# bin/import_dta_private.mex${MEXSUFFIX}: ${COMMON_OBJS} matlab/dta_to_matlab.o matlab/import_dta_private.cpp
# 	${MEX} ${MEXFLAGS} ${COMMON_INCLUDE_DIR} ${MATLAB_INCLUDE_DIR} $^ -output bin/import_dta_private

# matlab/dta_to_matlab.o: matlab/dta_to_matlab.cpp
# 	${CXX} ${CXXFLAGS} ${COMMON_INCLUDE_DIR} ${MATLAB_INCLUDE_DIR} -c $< -o $@

# yajl_build/%.o: yajl_build/%.c
# 	${CC} ${CFLAGS} -c $< -o $@

# common/%.o: common/%.cpp
# 	${CXX} ${CXXFLAGS} -c $< -o $@

# json/%.o: json/%.cpp
# 	${CXX} ${CXXFLAGS} ${COMMON_INCLUDE_DIR} ${YAJL_INCLUDE_DIR} -c $< -o $@

# clean:
# 	rm -rf *.o
# 	rm -rf common/*.o
# 	rm -rf yajl_build/*.o
# 	rm -rf json/*.o
# 	rm -rf matlab/*.o
# 	rm -rf *.a
# 	rm -f json/dta2json
# 	rm -f matlab/*.mex${MEXSUFFIX}