#!/bin/bash

##############################
## Welcome to SDC analysis! ##
##############################
# run with -h flag to see usage
# @author Atamert Rahma rahma@in.tum.de


usage() {
    echo "Usage: $0 [-b <path-to-build-dir>] [-n <number of processes>] [-x <field size in x>] [-y <field size in y>] [-t <simulation time>] [-d decomposition factor] [-r totalRuns]"
    1>&2;
    exit 1;
}

# get user input parameters for swe tests
while getopts b:n:x:y:t:d:r: flag
do
    case "${flag}" in
        b) b=${OPTARG};;
        n) np=${OPTARG};;
        x) x=${OPTARG};;
        y) y=${OPTARG};;
        t) t=${OPTARG};;
        d) d=${OPTARG};;
        r) totalRuns=${OPTARG};;
        *) usage;;
    esac
done
shift $((OPTIND-1))

if [ -z "${b}" ] || [ -z "${np}" ] || [ -z "${x}" ] || [ -z "${y}" ] || [ -z "${t}" ] || [ -z "${d}" ] || [ -z "${totalRuns}" ]; then
    usage
fi

if [ $TEAMS -ne 2 ]
then
    echo "set the environment variable TEAMS to 2 before running the tests!"
    exit 1
fi

#########################################
## RUN METHOD 1 FOR REFERENCE SOLUTION ##
#########################################
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

############################################
## RUN METHOD 3 AND 4 WHILE INJECTING SDC ##
############################################
# totalRuns = how many runs to make with each method, increase to get a more accurate percentage
numberOfInjectedSDC=$totalRuns
# number of corrected SDCs by method 3
method3_corrected=0
# number of failed runs by method 3 (exit code != 0)
method3_failed=0
# number of corrected SDCs by method 4
method4_corrected=0
# number of failed runs by method 4 (exit code != 0)
method4_failed=0

while [ $totalRuns -ne 0 ]
do
    #   --  METHOD 3 : softRes_admiss_useShared --  ##
    let np_3=$np
    outputPrefix_3='TEST_softRes_admiss_useShared_t'
    runPrefix_3='TEST_softRes_admiss_useShared'
    echo "mpirun --oversubscribe -np $np_3 $b/swe_softRes_admiss_useShared -x $x -y $y -t $t -o $runPrefix_3 -w -i 1 -d $d -f $(( $t / 2 ))"
    eval "mpirun --oversubscribe -np $np_3 $b/swe_softRes_admiss_useShared -x $x -y $y -t $t -o $runPrefix_3 -w -i 1 -d $d -f $(( $t / 2 ))"
    if [ $? -ne 0 ]
    then
        echo "ERROR: METHOD 3 : softRes_admiss_useShared has failed!"
        method3_failed=$(( $method3_failed + 1 ))
    fi

    #   --  METHOD 4 : softRes_admiss_redundant --  ##
    let np_4=2*$np
    outputPrefix_4='TEST_softRes_admiss_redundant_t'
    runPrefix_4='TEST_softRes_admiss_redundant'
    echo "mpirun --oversubscribe -np $np_4 $b/swe_softRes_admiss_redundant -x $x -y $y -t $t -o $runPrefix_4 -w -i 1 -d $d -f $(( $t / 2 ))"
    eval "mpirun --oversubscribe -np $np_4 $b/swe_softRes_admiss_redundant -x $x -y $y -t $t -o $runPrefix_4 -w -i 1 -d $d -f $(( $t / 2 ))"
    if [ $? -ne 0 ]
    then
        echo "ERROR: METHOD 4 : softRes_admiss_redundant has failed!"
        method4_failed=$(( $method4_failed + 1 ))
    fi

    ########################################################
    ## COMPARE THE OUTPUTS AGAINST THE REFERENCE SOLUTION ##
    ########################################################
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

    method3Correct=1
    method3_t0_correct=1
    method3_t1_correct=1
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
                eval "cmp $f1 $f2"
                exitCode=$?
                if [ $exitCode -eq 0 ]
                then
                    #echo "OUTPUTS EQUAL: $f1 & $f2";
                    echo "";
                elif [ $exitCode -eq 1 ]
                then
                    if [ $i -eq 0 ]
                    then
                        method3_t0_correct=0
                    elif [ $i -eq 1 ]
                    then
                        method3_t1_correct=0
                    fi
                    method3Correct=0
                    echo "method 3 could not detect this SDC.."
                else
                    echo "ERROR: FILE NOT FOUND: $f1 or $f2 doesn't exists!"
                    exit 1;
                fi
            done
        done
    done
    if [ $method3Correct -eq 1 ]
    then
        method3_corrected=$(( $method3_corrected + 1 ))
    else
        anyTeamCorrected=0
        if [ $method3_t0_correct -eq 1 ]
        then
            echo "METHOD 3 TEAM 0 HAS CORRECTED THE INJECTED SDC!"
            anyTeamCorrected=1
        fi
        if [ $method3_t1_correct -eq 1 ]
        then
            echo "METHOD 3 TEAM 1 HAS CORRECTED THE INJECTED SDC!"
            anyTeamCorrected=1
        fi
        if [ $anyTeamCorrected -eq 1 ]
        then
            method3_corrected=$(( $method3_corrected + 1 ))
        fi
    fi

    method4Correct=1
    method4_t0_correct=1
    method4_t1_correct=1
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
                eval "cmp $f1 $f2"
                exitCode=$?
                if [ $exitCode -eq 0 ]
                then
                    echo "";
                elif [ $exitCode -eq 1 ]
                then
                    if [ $i -eq 0 ]
                    then
                        method4_t0_correct=0
                    elif [ $i -eq 1 ]
                    then
                        method4_t1_correct=0
                    fi
                    method4Correct=0
                    echo "method 4 could not detect this SDC.."
                else
                    echo "ERROR: FILE NOT FOUND: $f1 or $f2 doesn't exists!"
                    exit 1;
                fi
            done
        done
    done

    if [ $method4Correct -eq 1 ]
    then
        method4_corrected=$(( $method4_corrected + 1 ))
    else
        anyTeamCorrected=0
        if [ $method4_t0_correct -eq 1 ]
        then
            echo "METHOD 4 TEAM 0 HAS CORRECTED THE INJECTED SDC!"
            anyTeamCorrected=1
        fi
        if [ $method4_t1_correct -eq 1 ]
        then
            echo "METHOD 4 TEAM 1 HAS CORRECTED THE INJECTED SDC!"
            anyTeamCorrected=1
        fi
        if [ $anyTeamCorrected -eq 1 ]
        then
            method4_corrected=$(( $method4_corrected + 1 ))
        fi
    fi
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
    totalRuns=$(( $totalRuns - 1 ))
    echo "METHOD 3: corrected $method3_corrected of $(( $numberOfInjectedSDC - $totalRuns )) SDCs"
    echo "METHOD 3: failed $method3_failed times (UNKNOWN ERROR)"
    echo "METHOD 4: corrected $method4_corrected of $(( $numberOfInjectedSDC - $totalRuns )) SDCs"
    echo "METHOD 4: failed $method4_failed times (UNKNOWN ERROR)"
    echo "$totalRuns injections to finish.."
done

# remove noRes output
for (( j=0; j<$ranksPerTeam_1; j++ ))
do
    let localBlockPositionX=$j/$blockCountY_1;
    let localBlockPositionY=$j%$blockCountY_1;
    echo "rm ${outputPrefix_1}0_${localBlockPositionX}_$localBlockPositionY.nc";
    eval "rm ${outputPrefix_1}0_${localBlockPositionX}_$localBlockPositionY.nc";
done

echo "------- Analysis FINNISHED -------"
echo "METHOD 3: corrected $method3_corrected of $numberOfInjectedSDC SDCs"
echo "METHOD 3: failed $method3_failed times (UNKNOWN ERROR)"
echo "METHOD 4: corrected $method4_corrected of $numberOfInjectedSDC SDCs"
echo "METHOD 4: failed $method4_failed times (UNKNOWN ERROR)"
echo "==> M3 has a correction rate of $method3_corrected/$numberOfInjectedSDC = $(echo "scale=4; $method3_corrected/$numberOfInjectedSDC" | bc)"
echo "==> M4 has a correction rate of $method4_corrected/$numberOfInjectedSDC = $(echo "scale=4; $method4_corrected/$numberOfInjectedSDC" | bc)"
echo "----------------------------------"

exit 0
