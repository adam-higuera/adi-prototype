#!/bin/bash
#
#PBS -l walltime=0:15:0
#PBS -l nodes=3:ppn=3
#PBS -q lazy
#PBS -N emsquare2d50_9
#PBS -o emsq2d25procs.out
#PBS -e emsq2d25procs.err
#PBS -V


module load intel/intel-13.0.0
module load openmpi/openmpi-1.6.4_intel-13.0.0_torque-4.2.3_ib
module load hpctoolkit/hpctoolkit-5.3.2_openmpi-1.7.0_gcc-4.8.0_torque-4.1.4
LD_LIBRARY_PATH=/projects/adhi1756/boost_1_53_0/stage/lib:$LD_LIBRARY_PATH

cd /projects/adhi1756/adi-prototype

PROG=${PROG_NAME##*/}
MEASUREMENTS_DIR=${RESULTS_DIR}/hpctoolkit-${PROG}-measurements-$PBS_JOBID

# hpcstruct -o ${RESULTS_DIR}/emsq2d-${PBS_JOBID}.hpcstruct ${PROG_NAME}
# mpiexec hpcrun -t -o ${MEASUREMENTS_DIR} ${PROG_NAME} --bynodes > ${RESULTS_DIR}/timing_s${DOMAIN_SIZE}timing_p$((${PROCS_PER_EDGE}*${PROCS_PER_EDGE}))_t${NUM_STEPS}
mpiexec ${PROG_NAME} --bynodes > ${RESULTS_DIR}/str${DOMAIN_SIZE}timing$((${PROCS_PER_EDGE}*${PROCS_PER_EDGE}))

# hpcprof-mpi -o ${RESULTS_DIR}/hpctoolkit-p${PROCS_PER_EDGE}-s${DOMAIN_SIZE}-${PBS_JOBID} -S ${RESULTS_DIR}/emsq2d-${PBS_JOBID}.hpcstruct ${MEASUREMENTS_DIR}

# rm -f ${RESULTS_DIR}/emsq2d-${PBS_JOBID}.hpcstruct
# rm -rf ${MEASUREMENTS_DIR}