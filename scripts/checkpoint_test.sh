#!/bin/bash
unset HOST
unset HOSTNAME

HOST=$(hostname)
echo $HOST

NODES=$(scontrol show hostnames ${SLURM_JOB_NODELIST})
NODES=(${NODES})
echo $HOSTNAME

APPLICATION="../build_noteam/swe-mpi"
MPI_PARAM=""
OUTPUT="log.txt"

PROCS=$SLURM_NTASKS
SIZE=${1:-4000}
MTBF=${2:-30}
FAILS=${3:-0}
HEARTBEAT=${4:-30}
RANDOM=0

export SPARES=$NUM_SPARES

echo "SIZE: $SIZE, SPARES: $SPARES, MTBF: $MTBF, FAILS: $FAILS"

START=$(date +"%s")
$APPLICATION -x $SIZE -y $SIZE -o $SCRATCH/output/test1 -b $SCRATCH/backup/test1 -i $HEARTBEAT &
sleep $MTBF

for i in $(seq 1 $FAILS); do
    pids=($(pgrep swe-mpi))
    num_nodes=${#NODES[@]}
    num_pids=${#pids[@]}

    if (($num_pids == 0)); then
        echo "This should not happen"
        break
    fi

    echo "Killing SWE with pid ${pids[0]}"
    kill -SIGKILL ${pids[0]}

    while true; do
        sleep 1
        pids=($(pgrep swe-mpi))
        num_pids=${#pids[@]}
        if (($num_pids == 0)); then
       	   break
   	fi
    done 

    $APPLICATION -x $SIZE -y $SIZE -o $SCRATCH/output/test1 -b $SCRATCH/backup/test1 -i $HEARTBEAT -r $SCRATCH/backup/test1 & 
    sleep $MTBF

done

echo "Waiting until SWE terminates"

while true; do
    sleep 1
    pids=($(pgrep swe-mpi))
    num_pids=${#pids[@]}
    if (($num_pids == 0)); then
        break
    fi
done

echo "SWE terminated"
END=$(date +"%s")
DURATION=$((END - START))

if [ "$SLURM_PROCID" = "0" ]; then
    echo "SIZE: $SIZE, SPARES: $NUM_SPARES, PROCS: $PROCS, MTBF: $MTBF, FAILS: $FAILS, DURATION: $DURATION, CP_INT: $HEARTBEAT" >>"teaMPI_log_cr.txt"
fi
