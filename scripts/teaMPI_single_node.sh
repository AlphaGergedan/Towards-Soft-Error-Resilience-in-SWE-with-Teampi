APPLICATION="../build/swe-mpi"
MPI_PARAM=""
OUTPUT="log.txt"

SIZE=4000
NUM_SPARES=4
MTBF=30
HEARTBEAT=5
FAILS=0

export SPARES=$NUM_SPARES

echo "$APPLICATION -x $SIZE -y $SIZE -o ../build/output/test1 -b ../build/backup/test1 -i $HEARTBEAT"

START=$(date +"%s")
$APPLICATION -x $SIZE -y $SIZE -o ../build/output/test1 -b ../build/backup/test1 -i $HEARTBEAT &
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
