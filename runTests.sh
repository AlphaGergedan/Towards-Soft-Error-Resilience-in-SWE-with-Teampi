#!/bin/bash

#############################
## Welcome to runTests.sh! ##
#############################
# run with -h flag to see usage
# @author Atamert Rahma rahma@in.tum.de


usage() {
    echo "Usage: $0 [-b <path-to-build-dir>] [-n <number of processes>] [-x <field size in x>] [-y <field size in y>] [-t <simulation time>] [-d decomposition factor]"
    1>&2;
    exit 1;
}

# get user input parameters for swe tests
while getopts b:n:x:y:t:d: flag
do
    case "${flag}" in
        b) b=${OPTARG};;
        n) np=${OPTARG};;
        x) x=${OPTARG};;
        y) y=${OPTARG};;
        t) t=${OPTARG};;
        d) d=${OPTARG};;
        *) usage;;
    esac
done
shift $((OPTIND-1))

if [ -z "${b}" ] || [ -z "${np}" ] || [ -z "${x}" ] || [ -z "${y}" ] || [ -z "${t}" ] || [ -z "${d}" ]; then
    usage
fi

if [ $TEAMS -ne 2 ]
then
    echo "set the environment variable TEAMS to 2 before running the tests!"
    exit 1
fi

####################################################################################
## RUN ALL THE METHODS ONCE WITHOUT INJECTING SDC | SEE IF THEY ARE EXITING WITH 0 ##
####################################################################################

#   --  METHOD 1 : noRes    --  ##
np_1=$(( np * d ))
echo "np is given as $np_1"
outputPrefix_1='TEST_noRes_t'
runPrefix_1='TEST_noRes'
echo "mpirun --oversubscribe -np $np_1 $b/swe_noRes -x $x -y $y -t $t -o $runPrefix_1 -w"
eval "mpirun --oversubscribe -np $np_1 $b/swe_noRes -x $x -y $y -t $t -o $runPrefix_1 -w"
if [ $? -eq 0 ]
then
    echo "-----------------------------------------------"
    echo "-- METHOD 1 : noRes is successfully finished --"
    echo "-----------------------------------------------"
else
    echo "ERROR: METHOD 1 : noRes has failed!"
    exit 1;
fi

ranksPerTeam_1=$np_1
blocksPerRank_1=1
totalBlocks_1=$ranksPerTeam_1
temp_1=$(echo "scale=4; sqrt($totalBlocks_1)" | bc)
blockCountY_1=${temp_1%%.*}
while [ $(( $totalBlocks_1 % $blockCountY_1 )) -ne 0 ]
do
  blockCountY_1=$(( $blockCountY_1 - 1 ))
done
let "blockCountX_1 = $totalBlocks_1 / $blockCountY_1"

for (( j=0; j<$ranksPerTeam_1; j++ ))
do
    let localBlockPositionX=$j/$blockCountY_1;
    let localBlockPositionY=$j%$blockCountY_1;
    echo "-> ${outputPrefix_1}0_${localBlockPositionX}_$localBlockPositionY.nc is written";
done

#   --  METHOD 2 : softRes_hashes --  ##
# run swe_softRes_hashes with -np 2*${np}*${d}
let np_2=2*$np*$d
outputPrefix_2='TEST_softRes_hashes_t'
runPrefix_2='TEST_softRes_hashes'
echo "mpirun --oversubscribe -np $np_2 $b/swe_softRes_hashes -x $x -y $y -t $t -o $runPrefix_2 -w -c 5"
eval "mpirun --oversubscribe -np $np_2 $b/swe_softRes_hashes -x $x -y $y -t $t -o $runPrefix_2 -w -c 5"
if [ $? -eq 0 ]
then
    echo "--------------------------------------------------------"
    echo "-- METHOD 2 : softRes_hashes is successfully finished --"
    echo "--------------------------------------------------------"
else
    echo "ERROR: METHOD 2 : softRes_hashes has failed!";
    exit 1;
fi

let ranksPerTeam_2=$np_2/$TEAMS
blocksPerRank_2=1
let totalBlocks_2=$blocksPerRank_2*$ranksPerTeam_2
temp_2=$(echo "scale=4; sqrt($totalBlocks_2)" | bc)
blockCountY_2=${temp_2%%.*}
while [ $(( $totalBlocks_2 % $blockCountY_2 )) -ne 0 ]
do
  blockCountY_2=$(( $blockCountY_2 - 1 ))
done
let "blockCountX_2 = $totalBlocks_2 / $blockCountY_2"

for (( i=0; i<$TEAMS; i++ ))
do
    for (( j=0; j<$ranksPerTeam_2; j++ ))
    do
        let localBlockPositionX=$j/$blockCountY_2;
        let localBlockPositionY=$j%$blockCountY_2;
        echo "-> ${outputPrefix_2}${i}_${localBlockPositionX}_$localBlockPositionY.nc is written";
    done
done

#   --  METHOD 3 : softRes_admiss_useShared --  ##
let np_3=$np
outputPrefix_3='TEST_softRes_admiss_useShared_t'
runPrefix_3='TEST_softRes_admiss_useShared'
echo "mpirun --oversubscribe -np $np_3 $b/swe_softRes_admiss_useShared -x $x -y $y -t $t -o $runPrefix_3 -w -i 1 -d $d"
eval "mpirun --oversubscribe -np $np_3 $b/swe_softRes_admiss_useShared -x $x -y $y -t $t -o $runPrefix_3 -w -i 1 -d $d"
if [ $? -eq 0 ]
then
    echo "------------------------------------------------------------------"
    echo "-- METHOD 3 : softRes_admiss_useShared is successfully finished --"
    echo "------------------------------------------------------------------"
else
    echo "ERROR: METHOD 3 : softRes_admiss_useShared has failed!"
    exit 1;
fi

let ranksPerTeam_3=$np_3/$TEAMS
let blocksPerRank_3=$TEAMS*$d
let totalBlocks_3=$blocksPerRank_3*$ranksPerTeam_3
temp_3=$(echo "scale=4; sqrt($totalBlocks_3)" | bc)
blockCountY_3=${temp_3%%.*}
while [ $(( $totalBlocks_3 % $blockCountY_3 )) -ne 0 ]
do
  blockCountY_3=$(( $blockCountY_3 - 1 ))
done
let "blockCountX_3 = $totalBlocks_3 / $blockCountY_3"

for (( i=0; i<$TEAMS; i++ ))
do
    for (( j=0; j<$ranksPerTeam_3; j++ ))
    do
        startPoint=$(( $j * $blocksPerRank_3 ))
        for (( currentBlockNr=$startPoint; currentBlockNr<$(( $startPoint + $blocksPerRank_3 )); currentBlockNr++ ))
        do
            let localBlockPositionX=$currentBlockNr/$blockCountY_3
            let localBlockPositionY=$currentBlockNr%$blockCountY_3;
            echo "-> ${outputPrefix_3}${i}_${localBlockPositionX}_$localBlockPositionY.nc is written";
        done
    done
done

#   --  METHOD 4 : softRes_admiss_redundant --  ##
let np_4=2*$np
outputPrefix_4='TEST_softRes_admiss_redundant_t'
runPrefix_4='TEST_softRes_admiss_redundant'
echo "mpirun --oversubscribe -np $np_4 $b/swe_softRes_admiss_redundant -x $x -y $y -t $t -o $runPrefix_4 -w -i 1 -d $d"
eval "mpirun --oversubscribe -np $np_4 $b/swe_softRes_admiss_redundant -x $x -y $y -t $t -o $runPrefix_4 -w -i 1 -d $d"
if [ $? -eq 0 ]
then
    echo "------------------------------------------------------------------"
    echo "-- METHOD 4 : softRes_admiss_redundant is successfully finished --"
    echo "------------------------------------------------------------------"
else
    echo "ERROR: METHOD 4 : softRes_admiss_redundant has failed!"
    exit 1;
fi

let ranksPerTeam_4=$np_4/$TEAMS
let blocksPerRank_4=$d
let totalBlocks_4=$blocksPerRank_4*$ranksPerTeam_4
temp_4=$(echo "scale=4; sqrt($totalBlocks_4)" | bc)
blockCountY_4=${temp_4%%.*}
while [ $(( $totalBlocks_4 % $blockCountY_4 )) -ne 0 ]
do
  blockCountY_4=$(( $blockCountY_4 - 1 ))
done
let "blockCountX_4 = $totalBlocks_4 / $blockCountY_4"

for (( i=0; i<$TEAMS; i++ ))
do
    for (( j=0; j<$ranksPerTeam_4; j++ ))
    do
        startPoint=$(( $j * $blocksPerRank_4 ))
        for (( currentBlockNr=$startPoint; currentBlockNr<$(( $startPoint + $blocksPerRank_4 )); currentBlockNr++ ))
        do
            let localBlockPositionX=$currentBlockNr/$blockCountY_4
            let localBlockPositionY=$currentBlockNr%$blockCountY_4;
            echo "-> ${outputPrefix_4}${i}_${localBlockPositionX}_$localBlockPositionY.nc is written";
        done
    done
done

###################################################
## RUNS FINISHED NOW WE WILL COMPARE THE RESULTS ##
###################################################

echo ""
echo " -- RUNS ARE FINISHED | NOW WE WILL COMPARE THE RESULTS -- "
echo ""

sleep 1

method2Correct=1

# compare method 2 with method 1
for (( i=0; i<$TEAMS; i++ ))
do
    for (( j=0; j<$ranksPerTeam_2; j++ ))
    do
        let localBlockPositionX=$j/$blockCountY_2;
        let localBlockPositionY=$j%$blockCountY_2;
        f1="${outputPrefix_1}0_${localBlockPositionX}_$localBlockPositionY.nc"
        f2="${outputPrefix_2}${i}_${localBlockPositionX}_$localBlockPositionY.nc"
        echo "--> cmp $f1 $f2"
        eval "cmp $f1 $f2"
        if [ $? -eq 0 ]
        then
            #echo "OUTPUTS EQUAL: $f1 & $f2";
            echo "";
        elif [ $? -eq 1 ]
        then
            echo "ERROR: OUTPUTS NOT EQUAL: $f1 & $f2 DIFFER!";
            method2Correct=0
        else
            echo "ERROR: FILE NOT FOUND: $f1 or $f2 doesn't exists!"
            exit 1;
        fi
    done
done

if [ $method2Correct -eq 0 ]
then
    echo "! WARNING: METHOD 2 PRODUCED A DIFFERENT OUTPUT !"
else
    echo "******************************************************************"
    echo "** METHOD 2 : method 2 produced the same output as the method 1 **"
    echo "******************************************************************"
fi

method3Correct=1
# compare method 3 with method 1
for (( i=0; i<$TEAMS; i++ ))
do
    for (( j=0; j<$ranksPerTeam_3; j++ ))
    do
        startPoint=$(( $j * $blocksPerRank_3 ))
        for (( currentBlockNr=$startPoint; currentBlockNr<$(( $startPoint + $blocksPerRank_3 )); currentBlockNr++ ))
        do
            let localBlockPositionX=$currentBlockNr/$blockCountY_3
            let localBlockPositionY=$currentBlockNr%$blockCountY_3;
            f1="${outputPrefix_1}0_${localBlockPositionX}_$localBlockPositionY.nc"
            f2="${outputPrefix_3}${i}_${localBlockPositionX}_$localBlockPositionY.nc"
            echo "--> cmp $f1 $f2"
            eval "cmp $f1 $f2"
            if [ $? -eq 0 ]
            then
                #echo "OUTPUTS EQUAL: $f1 & $f2";
                echo "";
            elif [ $? -eq 1 ]
            then
                echo "ERROR: OUTPUTS NOT EQUAL: $f1 & $f2 DIFFER!";
                method3Correct=0
            else
                echo "ERROR: FILE NOT FOUND: $f1 or $f2 doesn't exists!"
                exit 1;
            fi
        done
    done
done

if [ $method3Correct -eq 0 ]
then
    echo "! WARNING: METHOD 3 PRODUCED A DIFFERENT OUTPUT !"
else
    echo "******************************************************************"
    echo "** METHOD 3 : method 3 produced the same output as the method 1 **"
    echo "******************************************************************"
fi

method4Correct=1
# compare method 4 with method 1
for (( i=0; i<$TEAMS; i++ ))
do
    for (( j=0; j<$ranksPerTeam_4; j++ ))
    do
        startPoint=$(( $j * $blocksPerRank_4 ))
        for (( currentBlockNr=$startPoint; currentBlockNr<$(( $startPoint + $blocksPerRank_4 )); currentBlockNr++ ))
        do
            let localBlockPositionX=$currentBlockNr/$blockCountY_4
            let localBlockPositionY=$currentBlockNr%$blockCountY_4;
            f1="${outputPrefix_1}0_${localBlockPositionX}_$localBlockPositionY.nc"
            f2="${outputPrefix_4}${i}_${localBlockPositionX}_$localBlockPositionY.nc"
            echo "--> cmp $f1 $f2"
            eval "cmp $f1 $f2"
            if [ $? -eq 0 ]
            then
                #echo "OUTPUTS EQUAL: $f1 & $f2";
                echo "";
            elif [ $? -eq 1 ]
            then
                echo "ERROR: OUTPUTS NOT EQUAL: $f1 & $f2 DIFFER!";
                method4Correct=0
            else
                echo "ERROR: FILE NOT FOUND: $f1 or $f2 doesn't exists!"
                exit 1;
            fi
        done
    done
done

if [ $method4Correct -eq 0 ]
then
    echo "! WARNING: METHOD 4 PRODUCED A DIFFERENT OUTPUT !"
else
    echo "******************************************************************"
    echo "** METHOD 4 : method 4 produced the same output as the method 1 **"
    echo "******************************************************************"
fi

# remove noRes output
for (( j=0; j<$ranksPerTeam_1; j++ ))
do
    let localBlockPositionX=$j/$blockCountY_1;
    let localBlockPositionY=$j%$blockCountY_1;
    echo "rm ${outputPrefix_1}0_${localBlockPositionX}_$localBlockPositionY.nc";
    eval "rm ${outputPrefix_1}0_${localBlockPositionX}_$localBlockPositionY.nc";
done
# remove method 4 output
for (( i=0; i<$TEAMS; i++ ))
do
    for (( j=0; j<$ranksPerTeam_4; j++ ))
    do
        startPoint=$(( $j * $blocksPerRank_4 ))
        for (( currentBlockNr=$startPoint; currentBlockNr<$(( $startPoint + $blocksPerRank_4 )); currentBlockNr++ ))
        do
            let localBlockPositionX=$currentBlockNr/$blockCountY_4
            let localBlockPositionY=$currentBlockNr%$blockCountY_4;
            f="${outputPrefix_4}${i}_${localBlockPositionX}_$localBlockPositionY.nc"
            rm $f
            eval "rm BACKUP_${outputPrefix_4}${i}_${localBlockPositionX}_${localBlockPositionY}_metadata"
        done
    done
done
# remove method 3 output
for (( i=0; i<$TEAMS; i++ ))
do
    for (( j=0; j<$ranksPerTeam_3; j++ ))
    do
        startPoint=$(( $j * $blocksPerRank_3 ))
        for (( currentBlockNr=$startPoint; currentBlockNr<$(( $startPoint + $blocksPerRank_3 )); currentBlockNr++ ))
        do
            let localBlockPositionX=$currentBlockNr/$blockCountY_3
            let localBlockPositionY=$currentBlockNr%$blockCountY_3;
            f="${outputPrefix_3}${i}_${localBlockPositionX}_$localBlockPositionY.nc"
            rm $f
            eval "rm BACKUP_${outputPrefix_3}${i}_${localBlockPositionX}_${localBlockPositionY}_metadata"
        done
    done
done
# remove method 2 output
for (( i=0; i<$TEAMS; i++ ))
do
    for (( j=0; j<$ranksPerTeam_2; j++ ))
    do
        let localBlockPositionX=$j/$blockCountY_2;
        let localBlockPositionY=$j%$blockCountY_2;
        f="${outputPrefix_2}${i}_${localBlockPositionX}_$localBlockPositionY.nc"
        rm $f
    done
done

echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------"

if [ $(( $method2Correct + $method3Correct + $method4Correct )) -ne 3 ]
then
    echo "==> TESTS FAILED!!! ONE METHOD CREATED A DIFFERENT OUTPUT O.O"
    if [ $method2Correct -eq 0 ]
    then
        echo "  * method 2 : soft error detection with hashes -> created DIFFERENT outputs ==> FAILED"
    else
        echo "  * method 2 : soft error detection with hashes -> created the same outputs ==> PASSED"
    fi
    if [ $method3Correct -eq 0 ]
    then
        echo "  * method 3 : soft error resilience with admissibility checks and using task sharing -> created DIFFERENT outputs ==> FAILED"
    else
        echo "  * method 3 : soft error resilience with admissibility checks and using task sharing -> created the same outputs ==> PASSED"
    fi
    if [ $method4Correct -eq 0 ]
    then
        echo "  * method 4 : soft error resilience with admissibility checks and redundant computation -> created DIFFERENT outputs ==> FAILED"
    else
        echo "  * method 4 : soft error resilience with admissibility checks and redundant computation -> created the same outputs ==> PASSED"
    fi
    exit 1
else
    echo "==> TESTS ARE FINISHED! ALL OUTPUTS ARE EQUAL"
    exit 0;
fi
