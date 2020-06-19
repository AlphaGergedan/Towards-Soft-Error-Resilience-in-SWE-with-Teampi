#!bin/bash

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

SIZE=1000
NUM_SPARES=1
PROCS=$SLURM_NTASKS
MTBF=30
HEARTBEAT=5
FAILS=1
RANDOM=0

export SPARES=$NUM_SPARES
export OMP_NUM_THREADS=7

echo "$APPLICATION -x $SIZE -y $SIZE -o ../build/output/test1 -b ../build/backup/test1 -i $HEARTBEAT"

START=$(date +"%s")
$APPLICATION -x $SIZE -y $SIZE -o ../build/output/test1 -b ../build/backup/test1 -i $HEARTBEAT &
sleep 20

for i in $(seq 1 $FAILS); do
  pids=($(pgrep swe-mpi))
  num_pids=${#pids[@]}
  num_nodes=${#NODES[@]}
  if (($num_pids == 0)); then
    echo "This should not happen"
    break
  fi

  fail_proc=$(( ( RANDOM % $num_pids ) ))
  #fail_node=$(( (RANDOM % $num_nodes) ))
  fail_node=0

  echo "Killing proc $fail_proc of $num_pids still running procs on node ${NODES[$fail_node]}"
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

echo "SIZE: $SIZE, SPARES: $NUM_SPARES, PROCS: $PROCS, MTBF: $MTBF, FAILS: $FAILS, DURATION: $DURATION" >> "teaMPI_log.txt"
