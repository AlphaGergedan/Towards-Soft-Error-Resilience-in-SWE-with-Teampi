#!/bin/bash
APPLICATION="../build/swe-mpi"
MPI_PARAM="--mca mpi_ft_detector false"
OUTPUT="log.txt"


for((PROCS=2; PROCS<=8; PROCS=PROCS*2))
do
    for((i=60; i<=480; i=2*i))
    do
        for((SIZE=500; SIZE<=4000; SIZE=SIZE+500))
        do
            echo "/opt/bin/mpiexec $MPI_PARAM -np $PROCS $APPLICATION
            -x $SIZE -y $SIZE -o ../build/output/test1 -b ../build/backup/test1 -i $i"
            START=$(date +"%s")

            /opt/bin/mpiexec $MPI_PARAM -np $PROCS $APPLICATION -x $SIZE -y $SIZE -o ../build/output/test1 -b ../build/backup/test1 -i $i

            END=$(date +"%s")
            DURATION=$((END-START))

            FILENAME="${PROCS}_${i}_${OUTPUT}"
            echo "($SIZE, $DURATION)"
            echo "($SIZE, $DURATION)" >> $FILENAME

        done
    done
done


