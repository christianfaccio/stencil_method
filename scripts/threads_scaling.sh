#!/bin/bash

echo "Submitting threads scaling tests"

NODES=1
TASKS_PER_NODE=1
GRID_SIZE="-x 10000 -y 10000" 
N_STEPS="-n 500"
INJ_FREQ="-f 50"
N_SOURCES="-e 4"
E_SOURCE="-E 1.0"
OUTPUT="-o 0"
PERIODIC="-p 0"
VERBOSE="-v 0"

for threads in 1 2 4 8 16 32 56 84 112; do
    JOB_NAME="omp_scale_${threads}t"

    export SLURM_CPUS_PER_TASK=${threads}

    sbatch --nodes=${NODES} \
           --ntasks-per-node=${TASKS_PER_NODE} \
           --cpus-per-task=${threads} \
           --job-name=${JOB_NAME} \
	   --export=ALL,PROGRAM_ARGS="${GRID_SIZE} ${N_STEPS} ${INJ_FREQ} ${N_SOURCES} ${E_SOURCE} ${OUTPUT} ${PERIODIC} ${VERBOSE}",TEST_TYPE="omp" \
           scripts/go_dcgp.sbatch
done

echo "All OpenMp jobs submitted"
