#!/bin/bash
#
#PBS -l walltime=0:10:0
#PBS -l nodes=3:ppn=3
#PBS -q lazy
#PBS -N emsquare2d50_9
#PBS -o emsq2d25procs.out
#PBS -e emsq2d25procs.err
#PBS -V

cd /scr_verus/avh/Development/adi-prototype

# mpiexec hpcrun -t -o ${MEASUREMENTS_DIR} ${PROG_NAME} --bynodes > ${RESULTS_DIR}/str${DOMAIN_SIZE}timing$((${PROCS_PER_EDGE}*${PROCS_PER_EDGE}))
/usr/local/mpi/bin/mpiexec ${PROG_NAME} --bynodes > ${RESULTS_DIR}/timing_s${DOMAIN_SIZE}_p$((${PROCS_PER_EDGE}*${PROCS_PER_EDGE}))