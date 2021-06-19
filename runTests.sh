#!/bin/bash

#############################
## Welcome to runTests.sh! ##
#############################


usage() {
    echo "Usage: $0 [-b <path-to-build-dir>] [-n <number of processes>] [-x <field size in x>] [-y <field size in y>] [-t <simulation time>] [-d decomposition factor]"
    1>&2;
    exit 1;
}

removeAllOutputs() {
    method1_removeOutputs;
    method2_removeOutputs;
    method3_removeOutputs;
    method4_removeOutputs;
}
method1_removeOutputs() {
    for (( j=0; j<$ranksPerTeam_1; j++ ))
    do
        let localBlockPositionX=$j/$blockCountY_1;
        let localBlockPositionY=$j%$blockCountY_1;
        out="${outputPrefix_1}0_${localBlockPositionX}_$localBlockPositionY.nc"
        #eval "rm $out"
    done
}
method2_removeOutputs() {
    for (( i=0; i<$TEAMS; i++ ))
    do
        for (( j=0; j<$ranksPerTeam_2; j++ ))
        do
            let localBlockPositionX=$j/$blockCountY_2;
            let localBlockPositionY=$j%$blockCountY_2;
            out="${outputPrefix_2}${i}_${localBlockPositionX}_$localBlockPositionY.nc"
            #eval "rm $out"
        done
    done
}
method3_removeOutputs() {
    for (( i=0; i<$TEAMS; i++ ))
    do
        for (( j=0; j<$ranksPerTeam_3; j++ ))
        do
            let localBlockPositionX=$j/$blockCountY_3;
            let localBlockPositionY=$j%$blockCountY_3;
            prefix="${outputPrefix_3}${i}_${localBlockPositionX}_$localBlockPositionY"
            out="${prefix}.nc"
            backupOut="BACKUP_${prefix}_metadata"
            #eval "rm $out"
            #eval "rm $backupOut"
        done
    done

}
method4_removeOutputs() {
    for (( i=0; i<$TEAMS; i++ ))
    do
        for (( j=0; j<$ranksPerTeam_4; j++ ))
        do
            let localBlockPositionX=$j/$blockCountY_4;
            let localBlockPositionY=$j%$blockCountY_4;
            prefix="${outputPrefix_4}${i}_${localBlockPositionX}_$localBlockPositionY"
            out="${prefix}.nc"
            backupOut="BACKUP_${prefix}_metadata"
            #eval "rm $out"
            #eval "rm $backupOut"
        done
    done
}

# compare results of all the methods with the method 1 (no resilience)
compare() {
    # compare method 2 with method 1
    for (( i=0; i<$TEAMS; i++ ))
    do
        for (( j=0; j<$ranksPerTeam_2; j++ ))
        do
            let localBlockPositionX=$j/$blockCountY_2;
            let localBlockPositionY=$j%$blockCountY_2;
            f1="${outputPrefix_1}0_${localBlockPositionX}_$localBlockPositionY.nc"
            f2="${outputPrefix_2}${i}_${localBlockPositionX}_$localBlockPositionY.nc"
            echo "-> comparing the outputs: $f1 & $f2"
            eval "cmp $f1 $f2"
            if [ $? -eq 0 ]
            then
                echo "OUTPUTS EQUAL: $f1 & $f2";
            else
                echo "ERROR: OUTPUTS NOT EQUAL: $f1 & $f2 DIFFER!";
                exit 1;
            fi
        done
    done

    echo "******************************************************************"
    echo "** METHOD 2 : method 2 produced the same output as the method 1 **"
    echo "******************************************************************"
    #method2_removeOutputs

    # compare method 3 with method 1
    for (( i=0; i<$TEAMS; i++ ))
    do
        for (( j=0; j<$ranksPerTeam_3; j++ ))
        do
            let localBlockPositionX=$j/$blockCountY_3;
            let localBlockPositionY=$j%$blockCountY_3;
            f1="${outputPrefix_1}0_${localBlockPositionX}_$localBlockPositionY.nc"
            f2="${outputPrefix_3}${i}_${localBlockPositionX}_$localBlockPositionY.nc"
            echo "-> comparing the outputs: $f1 & $f2"
            eval "cmp $f1 $f2"
            if [ $? -eq 0 ]
            then
                echo "OUTPUTS EQUAL: $f1 & $f2";
            else
                echo "ERROR: OUTPUTS NOT EQUAL: $f1 & $f2 DIFFER!";
                exit 1;
            fi
        done
    done

    echo "******************************************************************"
    echo "** METHOD 3 : method 3 produced the same output as the method 1 **"
    echo "******************************************************************"
    #method3_removeOutputs

    # compare method 4 with method 1
    for (( i=0; i<$TEAMS; i++ ))
    do
        for (( j=0; j<$ranksPerTeam_4; j++ ))
        do
            let localBlockPositionX=$j/$blockCountY_4;
            let localBlockPositionY=$j%$blockCountY_4;
            f1="${outputPrefix_1}0_${localBlockPositionX}_$localBlockPositionY.nc"
            f2="${outputPrefix_4}${i}_${localBlockPositionX}_$localBlockPositionY.nc"
            echo "-> comparing the outputs: $f1 & $f2"
            eval "cmp $f1 $f2"
            if [ $? -eq 0 ]
            then
                echo "OUTPUTS EQUAL: $f1 & $f2";
            else
                echo "ERROR: OUTPUTS NOT EQUAL: $f1 & $f2 DIFFER!";
                exit 1;
            fi
        done
    done

    echo "******************************************************************"
    echo "** METHOD 4 : method 4 produced the same output as the method 1 **"
    echo "******************************************************************"
    #method4_removeOutputs
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

####################################################################################
## RUN ALL THE METHODS ONCE WITHOUT INJECTING SDC | SEE IF THEY ARE EXITING WITH 0 ##
####################################################################################

#   --  METHOD 1 : noRes    --  ##
np_1=$(( np * d ))
echo "np is given as $np_1"
outputPrefix_1='TEST_noRes_t'
runPrefix_1='TEST_noRes'
echo "mpirun --oversubscribe --display-map -np $np_1 $b/swe_noRes -x $x -y $y -t $t -o $runPrefix_1 -w"
eval "mpirun --oversubscribe --display-map -np $np_1 $b/swe_noRes -x $x -y $y -t $t -o $runPrefix_1 -w"
if [ $? -eq 0 ]
then
    echo "-----------------------------------------------"
    echo "-- METHOD 1 : noRes is successfully finished --"
    echo "-----------------------------------------------"
else
    echo "ERROR: METHOD 1 : noRes has failed!"
    exit 1;
fi

ranksPerTeam_1=$np
blocksPerRank_1=1
totalBlocks_1=$ranksPerTeam_1
temp_1=$(echo "scale=4; sqrt($totalBlocks_1)" | bc)
blockCountY_1=${temp_1%%.*}
while [ $(( totalBlocks_1 % blockCountY_1 )) -ne 0 ]
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
echo "mpirun --oversubscribe --display-map -np $np_2 $b/swe_softRes_hashes -x $x -y $y -t $t -o $runPrefix_2 -w -c 5"
eval "mpirun --oversubscribe --display-map -np $np_2 $b/swe_softRes_hashes -x $x -y $y -t $t -o $runPrefix_2 -w -c 5"
if [ $? -eq 0 ]
then
    echo "--------------------------------------------------------"
    echo "-- METHOD 2 : softRes_hashes is successfully finished --"
    echo "--------------------------------------------------------"
else
    echo "ERROR: METHOD 2 : softRes_hashes has failed!";
    exit 1;
fi

let ranksPerTeam_2=$np/$TEAMS
blocksPerRank_2=1
let totalBlocks_2=$blocksPerRank_2*$ranksPerTeam_2
temp_2=$(echo "scale=4; sqrt($totalBlocks_2)" | bc)
blockCountY_2=${temp_2%%.*}
while [ $(( totalBlocks_2 % blockCountY_2 )) -ne 0 ]
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
# run swe_softRes_admiss_useShared with -np ${np}   & compare results & delete results TODO fix the name of the executable after testing
let np_3=$np
outputPrefix_3='TEST_softRes_admiss_useShared_t'
runPrefix_3='TEST_softRes_admiss_useShared'
echo "mpirun --oversubscribe --display-map -np $np_3 $b/swe_softRes_admiss_useShared_v2 -x $x -y $y -t $t -o $runPrefix_3 -w -i 1 -d $d"
eval "mpirun --oversubscribe --display-map -np $np_3 $b/swe_softRes_admiss_useShared_v2 -x $x -y $y -t $t -o $runPrefix_3 -w -i 1 -d $d"
if [ $? -eq 0 ]
then
    echo "------------------------------------------------------------------"
    echo "-- METHOD 3 : softRes_admiss_useShared is successfully finished --"
    echo "------------------------------------------------------------------"
else
    echo "ERROR: METHOD 3 : softRes_admiss_useShared has failed!"
    exit 1;
fi

let ranksPerTeam_3=$np/$TEAMS
let blocksPerRank_3=$TEAMS*$d
let totalBlocks_3=$blocksPerRank_3*$ranksPerTeam_3
temp_3=$(echo "scale=4; sqrt($totalBlocks_3)" | bc)
blockCountY_3=${temp_3%%.*}
while [ $(( totalBlocks_3 % blockCountY_3 )) -ne 0 ]
do
  blockCountY_3=$(( $blockCountY_3 - 1 ))
done
let "blockCountX_3 = $totalBlocks_3 / $blockCountY_3"

for (( i=0; i<$TEAMS; i++ ))
do
    for (( j=0; j<$ranksPerTeam_3; j++ ))
    do
        let localBlockPositionX=$j/$blockCountY_3;
        let localBlockPositionY=$j%$blockCountY_3;
        echo "-> ${outputPrefix_3}${i}_${localBlockPositionX}_$localBlockPositionY.nc is written";
    done
done

#   --  METHOD 4 : softRes_admiss_redundant --  ##
# run swe_softRes_admiss_redundant with -np ${np}   & compare results & delete results
let np_4=$np
outputPrefix_4='TEST_softRes_admiss_redundant_t'
runPrefix_4='TEST_softRes_admiss_redundant'
echo "mpirun --oversubscribe --display-map -np $np_4 $b/swe_softRes_admiss_useShared_v2 -x $x -y $y -t $t -o $runPrefix_4 -w -i 1 -d $d"
eval "mpirun --oversubscribe --display-map -np $np_4 $b/swe_softRes_admiss_useShared_v2 -x $x -y $y -t $t -o $runPrefix_4 -w -i 1 -d $d"
if [ $? -eq 0 ]
then
    echo "------------------------------------------------------------------"
    echo "-- METHOD 4 : softRes_admiss_redundant is successfully finished --"
    echo "------------------------------------------------------------------"
else
    echo "ERROR: METHOD 4 : softRes_admiss_redundant has failed!"
    exit 1;
fi

let ranksPerTeam_4=$np/$TEAMS
let blocksPerRank_4=$d
let totalBlocks_4=$blocksPerRank_4*$ranksPerTeam_4
temp_4=$(echo "scale=4; sqrt($totalBlocks_4)" | bc)
blockCountY_4=${temp_4%%.*}
while [ $(( totalBlocks_4 % blockCountY_4 )) -ne 0 ]
do
  blockCountY_4=$(( $blockCountY_4 - 1 ))
done
let "blockCountX_4 = $totalBlocks_4 / $blockCountY_4"

for (( i=0; i<$TEAMS; i++ ))
do
    for (( j=0; j<$ranksPerTeam_4; j++ ))
    do
        let localBlockPositionX=$j/$blockCountY_4;
        let localBlockPositionY=$j%$blockCountY_4;
        echo "-> ${outputPrefix_4}${i}_${localBlockPositionX}_$localBlockPositionY.nc is written";
    done
done

###################################################
## RUNS FINISHED NOW WE WILL COMPARE THE RESULTS ##
###################################################
compare
echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
removeAllOutputs
exit 1;

# TODO also inject SDC and see if it is fixed

#########################################################################################
## RUN ALL THE METHODS AND INJECT SDC | WE WILL SEE THE PERCENTAGE OF CORRECTABLE SDCs ##
#########################################################################################

# how many runs to make with each method, increase to get a more accurate percentage
totalRuns=2

# number of detected SDCs by method 2
method2_detected=0
# number of failed runs by method 2 (exit code != 0)
method2_failed=0

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
    #   --  METHOD 2 : softRes_hashes --  ##
    # run swe_softRes_hashes with -np 2*${np}*${d}
    let np_2=2*$np*$d
    blockComparisonFailed_2=0
    outputPrefix_2='TEST_softRes_hashes_t'
    echo "mpirun --oversubscribe --display-map -np $np_2 $b/swe_softRes_hashes -x $x -y $y -t $t -o $outputPrefix_2 -w -c 5 -f $(( t / 2 ))"
    eval "mpirun --oversubscribe --display-map -np $np_2 $b/swe_softRes_hashes -x $x -y $y -t $t -o $outputPrefix_2 -w -c 5 -f $(( t / 2 ))"
    if [ $? -eq 0 ]
    then
        echo "---------------------------------------------------------------------------------------"
        echo "-- METHOD 2 : softRes_hashes finished with exit code 0. $totalRuns run remaining..  --"
        echo "---------------------------------------------------------------------------------------"
        # compare method 2 with method 1
        for (( i=0; i<$TEAMS; i++ ))
        do
            for (( j=0; j<$ranksPerTeam_2; j++ ))
            do
                let localBlockPositionX=$j/$blockCountY_2;
                let localBlockPositionY=$j%$blockCountY_2;
                f1="${outputPrefix_1}0_${localBlockPositionX}_$localBlockPositionY.nc"
                f2="${outputPrefix_2}${i}_${localBlockPositionX}_$localBlockPositionY.nc"
                echo "-> comparing the outputs: $f1 & $f2"
                eval "cmp $f1 $f2"
                if [ $? -eq 0 ]
                then
                    echo "OUTPUTS EQUAL: $f1 & $f2";
                    echo "BITFLIP DOESN'T HAVE AN EFFECT??";
                    echo "UNKNOWN ERROR: method 2 produced the same output under SDC";
                    exit 1;
                else
                    echo "ERROR: METHOD 2 COULDN'T DETECT AN SDC, POSSIBLY A HASH COLLISION!";
                    blockComparisonFailed_2=1;
                fi
            done
        done
    else
        echo "METHOD 2 DETECTED THE SDC!!";
        method2_detected=$(( $method2_detected + 1));
    fi
    #method2_removeOutputs
    if [ $blockComparisonFailed_2 -eq 1 ]
    then
        method2_failed=$(( $method2_failed + 1 ));
        blockComparisonFailed_2=0;
    fi

    #   --  METHOD 3 : softRes_admiss_useShared --  ##
    # run swe_softRes_admiss_useShared with -np ${np}   & compare results & delete results TODO fix the name of the executable after testing
    let np_3=$np
    blockComparisonFailed_3=0
    outputPrefix_3='TEST_softRes_admiss_useShared_t'
    echo "mpirun --oversubscribe --display-map -np $np_3 $b/swe_softRes_admiss_useShared_v2 -x $x -y $y -t $t -o $outputPrefix_3 -w -i 1 -d $d -f $(( t / 2 ))"
    eval "mpirun --oversubscribe --display-map -np $np_3 $b/swe_softRes_admiss_useShared_v2 -x $x -y $y -t $t -o $outputPrefix_3 -w -i 1 -d $d -f $(( t / 2 ))"
    if [ $? -eq 0 ]
    then
        echo "------------------------------------------------------------------------------------------------"
        echo "-- METHOD 3 : softRes_admiss_useShared finished with exit code 0. $totalRuns run remaining.. --"
        echo "------------------------------------------------------------------------------------------------"
        # compare method 3 with method 1
        for (( i=0; i<$TEAMS; i++ ))
        do
            for (( j=0; j<$ranksPerTeam_3; j++ ))
            do
                let localBlockPositionX=$j/$blockCountY_3;
                let localBlockPositionY=$j%$blockCountY_3;
                f1="${outputPrefix_1}0_${localBlockPositionX}_$localBlockPositionY.nc"
                f2="${outputPrefix_3}${i}_${localBlockPositionX}_$localBlockPositionY.nc"
                echo "-> comparing the outputs: $f1 & $f2"
                eval "cmp $f1 $f2"
                if [ $? -eq 0 ]
                then
                    echo "OUTPUTS EQUAL: $f1 & $f2";
                else
                    echo "ERROR: OUTPUTS NOT EQUAL: $f1 & $f2 DIFFER!";
                    blockComparisonFailed_3=1
                fi
            done
        done
    else
        echo "ERROR: METHOD 3 : softRes_admiss_useShared has failed!"
        method3_failed=$(( $method3_failed + 1))
    fi
    #method3_removeOutputs
    if [ $blockComparisonFailed_3 -eq 0 ]
    then
        method3_corrected=$(( $method3_corrected + 1 ));
    else
        echo "ERROR: METHOD 3 : COULDN'T CORRECT AN SDC!";
        blockComparisonFailed_3=0;
    fi

    #   --  METHOD 4 : softRes_admiss_redundant --  ##
    # run swe_softRes_admiss_redundant with -np ${np}   & compare results & delete results
    let np_4=$np
    blockComparisonFailed_4=0
    outputPrefix_4='TEST_softRes_admiss_redundant_t'
    echo "mpirun --oversubscribe --display-map -np $np_4 $b/swe_softRes_admiss_useShared_v2 -x $x -y $y -t $t -o $outputPrefix_4 -w -i 1 -d $d -f $(( t / 2 ))"
    eval "mpirun --oversubscribe --display-map -np $np_4 $b/swe_softRes_admiss_useShared_v2 -x $x -y $y -t $t -o $outputPrefix_4 -w -i 1 -d $d -f $(( t / 2 ))"
    if [ $? -eq 0 ]
    then
        echo "------------------------------------------------------------------------------------------------"
        echo "-- METHOD 4 : softRes_admiss_redundant finished with exit code 0. $totalRuns run remaining.. --"
        echo "------------------------------------------------------------------------------------------------"
        # compare method 4 with method 1
        for (( i=0; i<$TEAMS; i++ ))
        do
            for (( j=0; j<$ranksPerTeam_4; j++ ))
            do
                let localBlockPositionX=$j/$blockCountY_4;
                let localBlockPositionY=$j%$blockCountY_4;
                f1="${outputPrefix_1}0_${localBlockPositionX}_$localBlockPositionY.nc"
                f2="${outputPrefix_4}${i}_${localBlockPositionX}_$localBlockPositionY.nc"
                echo "-> comparing the outputs: $f1 & $f2"
                eval "cmp $f1 $f2"
                if [ $? -eq 0 ]
                then
                    echo "OUTPUTS EQUAL: $f1 & $f2";
                else
                    echo "ERROR: OUTPUTS NOT EQUAL: $f1 & $f2 DIFFER!";
                    blockComparisonFailed_4=1
                fi
            done
        done


    else
        echo "ERROR: METHOD 4 : softRes_admiss_redundant has failed!"
        method4_failed=$(( $method4_failed + 1 ))
    fi
    #method4_removeOutputs
    if [ $blockComparisonFailed_4 -eq 0 ]
    then
        method4_corrected=$(( $method4_corrected + 1 ));
    else
        echo "ERROR: METHOD 4 : COULDN'T CORRECT AN SDC!";
        blockComparisonFailed_4=0;
    fi

    totalRuns=$(( $totalRuns - 1 ))
done

echo "------- TESTS FINNISHED -------"
echo "METHOD 2: detected $method2_detected of 2 SDCs"
echo "METHOD 2: failed $method2_failed times (possibly a hash collision)"
echo "METHOD 3: corrected $method3_corrected of 2 SDCs"
echo "METHOD 3: failed $method3_failed times (UNKNOWN ERROR)"
echo "METHOD 4: corrected $method4_corrected of 2 SDCs"
echo "METHOD 4: failed $method4_failed times (UNKNOWN ERROR)"
echo "-------------------------------"

exit 0
