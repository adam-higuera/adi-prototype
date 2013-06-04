#!/bin/bash
#
#PBS -l walltime=0:10:0
#PBS -l nodes=3:ppn=3
#PBS -q lazy
#PBS -N emsquare2d50_9
#PBS -o emsq2d25procs.out
#PBS -e emsq2d25procs.err
#PBS -V

source /curc/tools/utils/dkinit

use .hpctoolkit_5.3.2_openmpi-1.6.3_intel-13.0.0_torque-2.5.11_ib
use .boost-1.50_openmpi-1.6_intel-12.1.4_torque-2.5.11_ib

cd /projects/adhi1756/adi-prototype

PROG=${PROG_NAME##*/}
MEASUREMENTS_DIR=${RESULTS_DIR}/hpctoolkit-${PROG}-measurements-$PBS_JOBID

hpcstruct -o ${RESULTS_DIR}/emsq2d-${PBS_JOBID}.hpcstruct ${PROG_NAME}
mpiexec hpcrun -t -o ${MEASUREMENTS_DIR} ${PROG_NAME} --bynodes > ${RESULTS_DIR}/timing_s${DOMAIN_SIZE}timing_p$((${PROCS_PER_EDGE}*${PROCS_PER_EDGE}))_t${NUM_STEPS}
# mpiexec ${PROG_NAME} --bynodes > ${RESULTS_DIR}/str${DOMAIN_SIZE}timing$((${PROCS_PER_EDGE}*${PROCS_PER_EDGE}))

hpcprof-mpi -o ${RESULTS_DIR}/hpctoolkit-p${PROCS_PER_EDGE}-s${DOMAIN_SIZE}-${PBS_JOBID} -S ${RESULTS_DIR}/emsq2d-${PBS_JOBID}.hpcstruct ${MEASUREMENTS_DIR}

rm -f ${RESULTS_DIR}/emsq2d-${PBS_JOBID}.hpcstruct
rm -rf ${MEASUREMENTS_DIR}