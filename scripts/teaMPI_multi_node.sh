#!/bin/bash

unset HOST
unset HOSTNAME

HOST=$(hostname)
echo $HOST 

NODES=$(scontrol show hostnames ${SLURM_JOB_NODELIST})
NODES=(${NODES})
echo $HOSTNAME

for ELEMENT in ${NODES[@]}; do
    echo "$HOST: knows $ELEMENT"
done

APPLICATION="../build/swe-mpi"
MPI_PARAM=""
OUTPUT="log.txt"

SIZE=${1:-4000}
NUM_SPARES=${2:-0}
PROCS=$SLURM_NTASKS
MTBF=${3:-30}
HEARTBEAT=5
FAILS=${4:-0}
RANDOM=0

export SPARES=$NUM_SPARES
export OMP_NUM_THREADS=28

echo "SIZE: $SIZE, SPARES: $SPARES, MTBF: $MTBF, FAILS: $FAILS"

START=$(date +"%s")
$APPLICATION -x $SIZE -y $SIZE -o ../build/output/test1 -b ../build/backup/test1 -i $HEARTBEAT &
sleep 20

for i in $(seq 1 $FAILS); do
  pids=($(pgrep swe-mpi))
  num_nodes=${#NODES[@]}
  if (($num_pids == 0)); then
    echo "This should not happen"
    break
  fi

  fail_node=$(( (RANDOM % $num_nodes) ))
  
  while (($fail_node == 0)); do 
    fail_node=$(( (RANDOM % $num_nodes) ))
  done

  echo "Killing SWE on node ${NODES[$fail_node]}"
  if [ "$HOST" = "${NODES[$fail_node]}" ]; then
    kill -SIGKILL ${pids[$fail_proc]}
  fi
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
DURATION=$((END-START))

if [ "$HOST" = "${NODES[0]}" ]; then
  echo "SIZE: $SIZE, SPARES: $NUM_SPARES, PROCS: $PROCS, MTBF: $MTBF, FAILS: $FAILS, DURATION: $DURATION" >> "teaMPI_log.txt"
fi
