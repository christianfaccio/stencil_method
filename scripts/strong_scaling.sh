#!/bin/bash

echo "Submitting Strong Scaling jobs..."

TASKS_PER_NODE=4
CPUS_PER_TASK=28
GRID_SIZE=10000 
N_STEPS=500
INJ_FREQ="-f 50"
N_SOURCES="-e 4"
E_SOURCE="-E 1.0"
OUTPUT="-o 0"
PERIODIC="-p 0"
VERBOSE="-v 0"

PARAMS="-x ${GRID_SIZE} -y ${GRID_SIZE} -n ${N_STEPS} ${INJ_FREQ} ${N_SOURCES} ${E_SOURCE} ${OUTPUT} ${PERIODIC} ${VERBOSE}"

for nodes in 1 2 4 8 16; do
    TOTAL_TASKS=$((nodes * TASKS_PER_NODE))
    JOB_NAME="strong_scale_${nodes}n_${TOTAL_TASKS}t"
    echo "Submitting job: ${JOB_NAME} with ${nodes} nodes and ${TOTAL_TASKS} tasks..."

    sbatch --nodes=${nodes} \
           --ntasks-per-node=${TASKS_PER_NODE} \
           --cpus-per-task=${CPUS_PER_TASK} \
           --job-name=${JOB_NAME} \
           --export=ALL,PROGRAM_ARGS="${PARAMS}",TEST_TYPE="strong" \
           scripts/go_dcgp.sbatch
done

echo "All Strong Scaling jobs submitted."
