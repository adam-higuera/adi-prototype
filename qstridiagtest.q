#!/bin/bash
#
#PBS -l walltime=4:0:0
#PBS -l nodes=5:ppn=5
#PBS -q lazy
#PBS -N tridiagtest
#PBS -o tridiagtest.out
#PBS -e tridiagtest.err
#PBS -V
#
cd /scr_verus/avh/Development/adi-prototype
export LD_LIBRARY_PATH=/scr_verus/avh/Development/boost_1_50_0/stage/lib:$LD_LIBRARY_PATH
#cd $PBS_O_WORKDIR
time /usr/local/mpi/bin/mpiexec tridiagtest