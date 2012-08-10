#!/bin/bash
#
#PBS -l walltime=4:0:0
#PBS -l nodes=12:ppn=8+1:ppn=4
#PBS -q immediate
#PBS -N emsquare2d
#PBS -o emsq2d25procs.out
#PBS -e emsq2d25procs.err
#PBS -V
#
cd /scr_verus/avh/Development/adi-prototype
export LD_LIBRARY_PATH=/scr_verus/avh/Development/boost_1_50_0/stage/lib:$LD_LIBRARY_PATH
#cd $PBS_O_WORKDIR
/usr/local/mpi/bin/mpiexec emsquare2d > timing100