#!/bin/bash
#SBATCH -J SWE
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=cm2
#SBATCH --partition=cm2_tiny
#SBATCH --qos=cm2_tiny
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=end
#SBATCH --mail-user=alexander.hoelzl@tum.de
#SBATCH --export=NONE
#SBATCH --time=03:00:00
  
module load slurm_setup
module unload intel-mpi
module load ulfm2/4.0.2u1-gcc8
module load netcdf-hdf5-all/4.7_hdf5-1.10-gcc8-serial   

unset HOST
unset HOSTNAME

HOST=$(hostname)
echo($HOST)

NODES=$(scontrol show hostnames ${SLURM_JOB_NODELIST})
NODES=(${NODES})
echo $HOSTNAME

APPLICATION="../build/swe-mpi"
MPI_PARAM=""
OUTPUT="log.txt"

SIZE=2000
NUM_SPARES=2
PROCS=$SLURM_NTASKS
MTBF=30
HEARTBEAT=5
FAILS=3

export SPARES=$NUM_SPARES

echo "mpiexec $MPI_PARAM -np $PROCS $APPLICATION -x $SIZE -y $SIZE -o ../build/output/test1 -b ../build/backup/test1 -i $HEARTBEAT"

START=$(date +"%s")
mpiexec $MPI_PARAM -np $PROCS $APPLICATION -x $SIZE -y $SIZE -o ../build/output/test1 -b ../build/backup/test1 -i $HEARTBEAT &
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
  fail_node=$(( (RANDOM % $num_nodes) ))

  if [ "$HOST" = "${NODES[fail_node]}"]; then
    echo "Killing proc $fail_proc of $num_pids still running procs on node ${NODES[fail_node]}"
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
