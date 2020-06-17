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
module load netcdf/4.7.2-gcc8-hdf5v1.10 

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

echo "/opt/bin/mpiexec $MPI_PARAM -np $PROCS $APPLICATION -x $SIZE -y $SIZE -o ../build/output/test1 -b ../build/backup/test1 -i $HEARTBEAT"

START=$(date +"%s")
/opt/bin/mpiexec $MPI_PARAM -np $PROCS $APPLICATION -x $SIZE -y $SIZE -o ../build/output/test1 -b ../build/backup/test1 -i $HEARTBEAT &
sleep 20

for i in $(seq 1 $FAILS); do
  pids=($(pgrep swe-mpi))
  num_pids=${#pids[@]}
  if (($num_pids == 0)); then
    echo "This should not happen"
    break
  fi
  fail_proc=$(( ( RANDOM % $num_pids ) ))
  echo "Killing proc $fail_proc of $num_pids still running procs"
  kill -SIGKILL ${pids[$fail_proc]}
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
