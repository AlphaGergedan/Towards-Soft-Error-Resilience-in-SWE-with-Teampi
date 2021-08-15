##----------------------------------------------------------------------##
##  SDC extraction script to decide if an injected SDC has resulted in  ##
##      - negligible error                                              ##
##      - correction                                                    ##
##      - DUE                                                           ##
##      - SDC                                                           ##
##  Note: this script needs to be run on the stdout output of the       ##
##        SDC injection script (runSDCAnalysis.sh), it should be run    ##
##        first.                                                        ##
##  Note: In this script M3 refers to the swe_softRes_admiss_useShared  ##
##        and M4 refers to the swe_softRes_admiss_redundant             ##
##        You can refer to README.md for more information about the     ##
##        methods.                                                      ##
##                                                                      ##
##  The script follows the following procedure:                         ##
##                                                                      ##
##                          SDC is injected                             ##
##                                |                                     ##
##                                |                                     ##
##                  SDC is reported in the application?                 ##
##                               /  \                                   ##
##                         (no) /    \ (yes)                            ##
##                             /      \                                 ##
##                outputs correct?    outputs correct?                  ##
##                       /  \               / \                         ##
##                 (no) /    \ (yes)  (no) /   \ (yes)                  ##
##                     /      \           /     \                       ##
##                    /        \         /       \                      ##
##              error in   NEGLIGIBLE  error in  CORRECTED              ##
##            application?            application?                      ##
##                / \                     / \                           ##
##          (no) /   \ (yes)        (no) /   \ (yes)                    ##
##              /     \                 /     \                         ##
##            SDC     DUE             SDC     DUE                       ##
##                                                                      ##
##  At the end of the script we get detection/correction rates of the   ##
##  methods against soft errors. In order to calcualate them we         ##
##  calculate totalInjections-NEGLIGIBLE and divide it with the number  ##
##  of total detections/corrections.                                    ##
##                                                                      ##
##  Usage: python extractSDC_outcomeRate.py <SDC-injection-output>   ##
##          <build-dir> <number-MPI-processes> <size-x> <size-y>        ##
##          <simulation-duration> <decomposition-factor>                ##
##  Author: Atamert Rahma | rahma@in.tum.de                             ##
##----------------------------------------------------------------------##

import sys

fileName = sys.argv[1]
given_build_dir = sys.argv[2]
given_np = sys.argv[3]
given_sizeX = sys.argv[4]
given_sizeY = sys.argv[5]
given_simTime = sys.argv[6]
given_decompfactor = sys.argv[7]

#f = open("TEAMS_SDC_injectionAnalysis_random_bitflip_n2_xy200_t100_d1_r1000.txt", 'r')
f = open (str(fileName), 'r')
call_useShared = "mpirun --oversubscribe -np " + given_np + " " + str(given_build_dir) + "/swe_softRes_admiss_useShared -x " + given_sizeX + " -y " + given_sizeY + " -t " + given_simTime + " -o TEST_softRes_admiss_useShared -w -i 1 -d " + given_decompfactor + " -f " + str(int (int(given_simTime)/2))
call_redundant = "mpirun --oversubscribe -np " + str(int(given_np)*2) + " " + str(given_build_dir) + "/swe_softRes_admiss_redundant -x " + given_sizeX + " -y " + given_sizeY + " -t " + given_simTime + " -o TEST_softRes_admiss_redundant -w -i 1 -d " + given_decompfactor + " -f " + str(int(int(given_simTime) / 2))

totalInjections = 0

M3_SDC_vanished = 0 #SDC was very small
M3_SDC_detected = 0 #error at the same time
M3_SDC_corrected = 0
M3_SDC = 0

M4_SDC_vanished = 0 #SDC was very small
M4_SDC_detected = 0 #error at the same time
M4_SDC_corrected = 0
M4_detected_otherTeamAlive = 0 #M4 has one healthy team! and it detected SDC!
M4_SDC = 0

block = ""
for line in f:
    # find the mpirun call to useShared
    if call_useShared in line:
        block += line
        indicator = "injections to finish.."
        # add lines until injection is finished
        for line2 in f:
            block += line2
            if indicator in line2:
                break

        blockLines = block.split('\n')
        # from mpi useShared call to redundant
        block_useShared = []
        for blockLine in blockLines:
            if call_redundant != blockLine:
                block_useShared.append(blockLine)
            else:
                break
        # from mpi redundant to end
        reached = False
        block_redundant = []
        for blockLine in blockLines:
            if reached:
                block_redundant.append(blockLine)
                continue
            else:
                if blockLine == call_redundant:
                    block_redundant.append(blockLine)
                    reached = True

        #M3 extraction
        M3_seems_error = False
        M3_reported_SDC = False
        M3error = "Errorhandler"
        M3_SDC_report = "SDC reported"
        for m3Line in block_useShared:
            if M3error in m3Line:
                M3_seems_error = True
            if M3_SDC_report in m3Line:
                M3_reported_SDC = True

        #M4 extraction
        M4_seems_error = False
        M4_reported_SDC = False
        M4_oneTeamAlive = False
        M4error = "Errorhandler"
        M4_SDC_report = "SDC reported"
        M4_teamCorrected_string = 'Redundant TEAM 1 HAS A CORRECT SOLUTION!'
        M3_seems_corrected = True
        M4_seems_corrected = True
        M3couldNotDetect = "Sharing could not detect this SDC.."
        M4couldNotDetect = "Redundant could not detect this SDC.."
        for m4Line in block_redundant:
            if M3couldNotDetect in m4Line:
                M3_seems_corrected = False
            if M4couldNotDetect in m4Line:
                M4_seems_corrected = False
            if M4error in m4Line:
                M4_seems_error = True
            if M4_SDC_report in m4Line:
                M4_reported_SDC = True
            if M4_teamCorrected_string in m4Line:
                M4_oneTeamAlive = True

        #M3
        if M3_reported_SDC:
            if M3_seems_corrected:
                M3_SDC_corrected += 1
            else:
                if M3_seems_error:
                    M3_SDC_detected += 1
        else:
            if M3_seems_corrected:
                M3_SDC_vanished += 1
            else:
                if M3_seems_error:
                    M3_SDC_detected += 1
                else:
                    M3_SDC += 1

        #M4
        if M4_reported_SDC:
            if M4_seems_corrected:
                M4_SDC_corrected += 1
            else:
                if M4_seems_error:
                    M4_SDC_detected += 1
                else:
                    if M4_oneTeamAlive:
                        M4_SDC_corrected += 1
                        M4_detected_otherTeamAlive += 1
                    else:
                        M4_SDC += 1
        else:
            if M4_seems_corrected:
                M4_SDC_vanished += 1
            else:
                if M4_seems_error:
                    M4_SDC_detected += 1
                else:
                    M4_SDC += 1
        block = ""
        totalInjections = totalInjections+1



nonVanishedInjections_M3 = totalInjections - M3_SDC_vanished
nonVanishedInjections_M4 = totalInjections - M4_SDC_vanished

perc_M3_vanished = round((M3_SDC_vanished / totalInjections), 4)
perc_M4_vanished = round((M4_SDC_vanished / totalInjections), 4)
perc_M3_detected = round((M3_SDC_detected / nonVanishedInjections_M3), 4)
perc_M4_detected = round((M4_SDC_detected / nonVanishedInjections_M4), 4)
perc_M3_corrected = round((M3_SDC_corrected / nonVanishedInjections_M3), 4)
perc_M4_corrected = round((M4_SDC_corrected / nonVanishedInjections_M4), 4)
perc_M4_oneTeamCorrected = round((M4_detected_otherTeamAlive / nonVanishedInjections_M4), 4)
perc_M3_SDC = round((M3_SDC / nonVanishedInjections_M3), 4)
perc_M4_SDC = round((M4_SDC / nonVanishedInjections_M4), 4)

if nonVanishedInjections_M3 != (M3_SDC + M3_SDC_detected + M3_SDC_corrected):
    print("WARNING: nonVanishedInjections_M3 = " + str(nonVanishedInjections_M3) + " != " + str(M3_SDC + M3_SDC_detected + M3_SDC_corrected) + " (M3_SDC + M3_SDC_detected + M3_SDC_corrected)")
if nonVanishedInjections_M4 != (M4_SDC + M4_SDC_detected + M4_SDC_corrected):
    print("WARNING: nonVanishedInjections_M4 = " + str(nonVanishedInjections_M4) + " != " + str(M4_SDC + M4_SDC_detected + M4_SDC_corrected) + " (M4_SDC + M4_SDC_detected + M4_SDC_corrected)")

print("**** Sharing **** " + str(totalInjections) + " total injections")
print(" -> SDC is VANISHED \t: " + str(M3_SDC_vanished) + "\t--> " + str(perc_M3_vanished) + "\tof total")
print(" -> SDC is DETECTED (ERROR): " + str(M3_SDC_detected) + "\t--> " + str(perc_M3_detected) + "\tof nonVanished")
print(" -> SDC is CORRECTED \t: " + str(M3_SDC_corrected) + "\t--> " + str(perc_M3_corrected) + "\tof nonVanished")
print(" -> SDC \t\t: " + str(M3_SDC) + "\t--> " + str(perc_M3_SDC) + "\tof nonVanished")
print(" -> correction rate = SDC CORRECTED / nonVanished  = " + str((M3_SDC_corrected) / nonVanishedInjections_M3))
print("**** Redundant **** " + str(totalInjections) + " total injections")
print(" -> SDC is VANISHED \t: " + str(M4_SDC_vanished) + "\t--> " + str(perc_M4_vanished) + "\tof total")
print(" -> SDC is DETECTED (ERROR): " + str(M4_SDC_detected) + "\t--> " + str(perc_M4_detected) + "\tof nonVanished")
print((" -> SDC is CORRECTED \t: " + str(M4_SDC_corrected)) + "\t--> " + str(perc_M4_corrected) + "\tof nonVanished " + (" | only 1 team corrected: " + str(M4_detected_otherTeamAlive) + "\t--> " + str(perc_M4_oneTeamCorrected) + "\tof nonVanished"))
print(" -> SDC \t\t: " + str(M4_SDC) + "\t--> " + str(perc_M4_SDC) + "\tof nonVanished")
print(" -> correction rate = SDC CORRECTED / nonVanished  = " + str((M4_SDC_corrected) / (M4_SDC + M4_SDC_detected + M4_SDC_corrected)))
print("************")
