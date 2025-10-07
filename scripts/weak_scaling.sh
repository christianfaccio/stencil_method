#!/bin/bash

echo "Submitting Weak Scaling jobs..."

TASKS_PER_NODE=8 # MPI tasks
CPUS_PER_TASK=14 # OMP threads
LOCAL_SIZE=4000 # single node grid size
N_STEPS=500
INJ_FREQ="-f 50"
N_SOURCES="-e 4"
E_SOURCE="-E 1.0"
OUTPUT="-o 0"
PERIODIC="-p 0"
VERBOSE="-v 0"

for nodes in 1 2 4 8 16; do
    TOTAL_TASKS=$((nodes * TASKS_PER_NODE))

    case ${TOTAL_TASKS} in
	8) PX=4; PY=2; ;;
        16) PX=4; PY=4; ;;
        32) PX=8; PY=4; ;;
        64) PX=8; PY=8; ;;
	128) PX=16; PY=8; ;;
        *) echo "Invalid number of tasks: ${TOTAL_TASKS}"; exit 1; ;;
    esac

    GRID_SIZE_X=$((PX * LOCAL_SIZE))
    GRID_SIZE_Y=$((PY * LOCAL_SIZE))
    PARAMS="-x ${GRID_SIZE_X} -y ${GRID_SIZE_Y} -n ${N_STEPS} ${INJ_FREQ} ${N_SOURCES} ${E_SOURCE} ${OUTPUT} ${PERIODIC} ${VERBOSE}"

    JOB_NAME="weak_scale_${nodes}n_${TOTAL_TASKS}t"
    echo "Submitting job: ${JOB_NAME} with ${nodes} nodes and ${TOTAL_TASKS} tasks..."

    sbatch --nodes=${nodes} \
           --ntasks-per-node=${TASKS_PER_NODE} \
           --cpus-per-task=${CPUS_PER_TASK} \
           --job-name=${JOB_NAME} \
           --export=ALL,PROGRAM_ARGS="${PARAMS}",TEST_TYPE="weak" \
           scripts/go_dcgp.sbatch
done

echo "All Weak Scaling jobs submitted."
