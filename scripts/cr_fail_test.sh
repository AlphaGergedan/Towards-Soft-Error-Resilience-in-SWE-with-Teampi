 #!/bin/bash
SWE_PATH="../build/swe-mpi"
MPI_PARAM="--mca mpi_ft_detector_thread true"
OUPUT="log.txt"
PROCS="4"
SIZE=3000

for((i=60; i<=480; i=2*i))
do
START=$(date +"%s")

/opt/bin/mpiexec $MPI_PARAM -np $PROCS $SWE_PATH -x $SIZE -y $SIZE -o ../build/output/test1 -b ../build/backup/test1 -i $i

END=$(date +"%s")
DURATION=$((END-START))

echo "($i, $DURATION)" >> "${PROCS}_${SIZE}_${OUPUT}"

done


 
