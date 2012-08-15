cmake_minimum_required(VERSION 2.6)
project(adi-prototype)

set(Boost_ADDITIONAL_VERSIONS "1.50" "1.50.0")

find_package(Boost COMPONENTS mpi program_options REQUIRED)
find_package(MPI REQUIRED)
find_package(LAPACK REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_DIRS})
include_directories(${LAPACK_INCLUDE_DIRS})

link_directories(${Boost_LIBRARY_DIRS})
link_directories(${MPI_LIBRARY_DIRS})
link_directories(${LAPACK_LIBRARY_DIRS})

add_executable(emsquare2d main.cpp MatrixBlock.cpp EMSquare2D.cpp)
add_executable(emsquare2d-strong strong_scaling.cpp MatrixBlock.cpp EMSquare2D.cpp)
add_executable(tridiagtest tridiagtest.cpp MatrixBlock.cpp)
add_executable(multitridiagtest multiple_tridiagtest.cpp MatrixBlock.cpp)

target_link_libraries(emsquare2d ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(emsquare2d-strong ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(tridiagtest ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(multitridiagtest ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${LAPACK_LIBRARIES})
