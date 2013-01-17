#!/bin/bash
#
#PBS -l walltime=0:5:0
#PBS -l nodes=3:ppn=3
#PBS -q lazy
#PBS -N emsquare2d50_9
#PBS -o emsq2d25procs.out
#PBS -e emsq2d25procs.err
#PBS -V

source /curc/tools/utils/dkinit

use .boost-1.50_openmpi-1.6_intel-12.1.4_torque-2.5.11_ib
use .hpctoolkit_5.2.1_openmpi-1.6_gcc-4.7.1_torque-2.5.11_ib

cd /projects/adhi1756/adi-prototype

mpiexec ${PROG_NAME} --bynodes > ${RESULTS_DIR}/str${DOMAIN_SIZE}timing$((${PROCS_PER_EDGE}*${PROCS_PER_EDGE}))
