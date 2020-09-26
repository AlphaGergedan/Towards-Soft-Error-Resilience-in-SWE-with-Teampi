#!/bin/bash
APPLICATION="../build_noteam/swe-mpi"

PROCS=$SLURM_NTASKS
SIZE=${1:-4000}
WAIT_TIME=${2:-30}
FAIL=${3:-0}
HEARTBEAT=${4:-30}
INIT_RUN=${5:-1}
RANDOM=0

echo "SIZE: $SIZE, WAIT_TIME: $WAIT_TIME, FAIL: $FAIL, HEARBEAT: $HEARTBEAT; INIT_RUN: $INIT_RUN"

#Start SWE either as fresh run or by loading from checkpoint
if (( $INIT_RUN )); then
    $APPLICATION -x $SIZE -y $SIZE -o $SCRATCH/output/test1 -b $SCRATCH/backup/test1 -i $HEARTBEAT &
    sleep $WAIT_TIME
else
    $APPLICATION -x $SIZE -y $SIZE -o $SCRATCH/output/test1 -b $SCRATCH/backup/test1 -i $HEARTBEAT -r $SCRATCH/backup/test1 &
    sleep $WAIT_TIME
fi

#If necessary inject failure
if (( $FAIL )); then

    pids=($(pgrep swe-mpi))
    num_nodes=${#NODES[@]}
    num_pids=${#pids[@]}

    if (( $num_pids == 0 )); then
        echo "This should not happen"
        break
    fi

    echo "Killing SWE with pid ${pids[0]}"
    kill -SIGKILL ${pids[0]}

fi

while true; do
   sleep 1
   pids=($(pgrep swe-mpi))
   num_pids=${#pids[@]}
   if (($num_pids == 0)); then
    break
   fi
done

echo "SWE terminated"