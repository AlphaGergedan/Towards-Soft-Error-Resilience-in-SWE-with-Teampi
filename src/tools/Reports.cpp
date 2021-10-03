/**
 * @file src/tools/Reports.cpp
 *
 * @brief Implementation of src/tools/Reports.hpp
 *
 * @author Atamert Rahma rahma@in.tum.de
 */

#include "tools/Reports.hpp"

tools::Reports::Reports(int numTeams, int myTeam, int myRankInTeam, unsigned int decompFactor, unsigned int blocksPerRank,
              std::vector<std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>> &simulationBlocks,
              std::vector<int> &myBlockOrder, unsigned char *primaryBlocksCorrupted,
              unsigned char *receivedBlocksCorrupted, unsigned char *SDC, unsigned char *SDC_inReplica,
              unsigned char *replicaCorrupted, unsigned char *receivedBlockReports,
              MPI_Comm interTeamComm) : simulationBlocks(simulationBlocks), myBlockOrder(myBlockOrder) {
    tools::Reports::numTeams = numTeams;
    tools::Reports::myTeam = myTeam;
    tools::Reports::myRankInTeam = myRankInTeam;
    tools::Reports::decompFactor = decompFactor;
    tools::Reports::blocksPerRank = blocksPerRank;
    tools::Reports::simulationBlocks = simulationBlocks;
    tools::Reports::myBlockOrder = myBlockOrder;
    tools::Reports::primaryBlocksCorrupted = primaryBlocksCorrupted;
    tools::Reports::receivedBlocksCorrupted = receivedBlocksCorrupted;
    tools::Reports::SDC = SDC;
    tools::Reports::SDC_inReplica = SDC_inReplica;
    tools::Reports::replicaCorrupted = replicaCorrupted;
    tools::Reports::receivedBlockReports = receivedBlockReports;
    tools::Reports::interTeamComm = interTeamComm;
}

tools::Reports::Reports(int numTeams, int myTeam, int myRankInTeam, unsigned int decompFactor, unsigned int blocksPerRank,
        std::vector<std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>> &simulationBlocks,
        unsigned char *primaryBlocksCorrupted, unsigned char *SDC,
        unsigned char *SDC_inReplica, unsigned char *replicaCorrupted, MPI_Comm interTeamComm)
    : simulationBlocks(simulationBlocks), myBlockOrder(myBlockOrder) {
    tools::Reports::numTeams = numTeams;
    tools::Reports::myTeam = myTeam;
    tools::Reports::myRankInTeam = myRankInTeam;
    tools::Reports::decompFactor = decompFactor;
    tools::Reports::blocksPerRank = blocksPerRank;
    tools::Reports::simulationBlocks = simulationBlocks;
    tools::Reports::primaryBlocksCorrupted = primaryBlocksCorrupted;
    tools::Reports::SDC = SDC;
    tools::Reports::SDC_inReplica = SDC_inReplica;
    tools::Reports::replicaCorrupted = replicaCorrupted;
    tools::Reports::interTeamComm = interTeamComm;
}

void tools::Reports::reportSDC() {
    /* report to all replicas */
    for (int destTeam = 0; destTeam < numTeams; destTeam++) {
        if (destTeam != myTeam)
            MPI_Send(SDC, 1, MPI_BYTE, destTeam, MPI_TAG_REPORT_PRIMARY_BLOCK, interTeamComm);
    }
}

int tools::Reports::getReloadReplica() {
    int reloadReplica;
    /* Receive the reload replica : tag 23 for reload replica */
    MPI_Recv(&reloadReplica, 1, MPI_INT, MPI_ANY_SOURCE,
             MPI_TAG_RECEIVE_RELOAD_REPLICA, interTeamComm, MPI_STATUS_IGNORE);
    return reloadReplica;
}

void tools::Reports::reportPrimaryBlocks(int reloadReplica) {
    /* send reports for all my primary blocks */
    for (unsigned int i = 0; i < decompFactor; i++)
        MPI_Send(primaryBlocksCorrupted+i, 1, MPI_BYTE, reloadReplica,
                 MPI_TAG_REPORT_PRIMARY_BLOCK, interTeamComm);
}

void tools::Reports::recoverCorruptedPrimaryBlocks(int reloadReplica, float t) {
    /* for each corrupted block receive b,h,hv,hu */
    for (unsigned int i = 0; i < decompFactor; i++) {
        if (primaryBlocksCorrupted[i] == 1) {
            /* Get a reference to the current corrupted block */
            const int& currentBlockNr = myBlockOrder[i];
            auto& currentCorruptedBlock = *simulationBlocks[currentBlockNr];
            /* Size of the arrays b,h,hv,hu */
            const int dataArraySize = currentCorruptedBlock.dataArraySize;
            std::cout << "T" << myTeam << "R" << myRankInTeam
                      << " : receiving b,h,hv,hu for my primaryBlock" << currentBlockNr
                      << ", block nr. is " << currentBlockNr
                      << std::endl;
            // primaryBlock recovery by tag 24 : receive order is important! --> b,h,hv,hu
            MPI_Recv(currentCorruptedBlock.b.getRawPointer(),
                     dataArraySize, MPI_FLOAT, reloadReplica,
                     MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                     interTeamComm, MPI_STATUS_IGNORE);
            MPI_Recv(currentCorruptedBlock.h.getRawPointer(),
                     dataArraySize, MPI_FLOAT, reloadReplica,
                     MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                     interTeamComm, MPI_STATUS_IGNORE);
            MPI_Recv(currentCorruptedBlock.hv.getRawPointer(),
                     dataArraySize, MPI_FLOAT, reloadReplica,
                     MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                     interTeamComm, MPI_STATUS_IGNORE);
            MPI_Recv(currentCorruptedBlock.hu.getRawPointer(),
                     dataArraySize, MPI_FLOAT, reloadReplica,
                     MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                     interTeamComm, MPI_STATUS_IGNORE);
            std::cout << "T" << myTeam << "R" << myRankInTeam
                      << " : b,h,hv,hu for my primaryBlock" << currentBlockNr
                      << " are received! Thanks replica " << reloadReplica
                      << std::endl;
            /* compute and validate again */
            currentCorruptedBlock.computeNumericalFluxes();
            bool admissible = currentCorruptedBlock.validateAdmissibility(t);
            if (!admissible)
                assert(false); // TODO we must abort right ?
            /* problem solved for this corrupted block */
            else {
                primaryBlocksCorrupted[i] = 0;
                std::cout << "T" << myTeam << "R" << myRankInTeam
                          << " : problem solved for my primaryBlock" << currentBlockNr
                          << ", block nr. is " << currentBlockNr
                          << std::endl;
                //MPI_Send(primaryBlocksCorrupted+i, 1, MPI_BYTE, reloadReplica, 100, interTeamComm); TODO we assume we solved the SDC
            }
        }
    }
}

void tools::Reports::recoverCorruptedPrimaryBlocks_redundant(int reloadReplica, float t) {
    /* for each corrupted block receive b,h,hv,hu */
    for (unsigned int i = 0; i < decompFactor; i++) {
        if (primaryBlocksCorrupted[i] == 1) {
            /* Get a reference to the current corrupted block */
            auto& currentCorruptedBlock = *simulationBlocks[i];
            /* Size of the arrays b,h,hv,hu */
            const int dataArraySize = currentCorruptedBlock.dataArraySize;
            std::cout << "T" << myTeam << "R" << myRankInTeam
                    << " : receiving b,h,hv,hu for my block " << i
                    << std::endl;
            // primaryBlock recovery by tag 24 : receive order is important! --> b,h,hv,hu
            MPI_Recv(currentCorruptedBlock.b.getRawPointer(),
                    dataArraySize, MPI_FLOAT, reloadReplica,
                    MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                    interTeamComm, MPI_STATUS_IGNORE);
            MPI_Recv(currentCorruptedBlock.h.getRawPointer(),
                    dataArraySize, MPI_FLOAT, reloadReplica,
                    MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                    interTeamComm, MPI_STATUS_IGNORE);
            MPI_Recv(currentCorruptedBlock.hv.getRawPointer(),
                    dataArraySize, MPI_FLOAT, reloadReplica,
                    MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                    interTeamComm, MPI_STATUS_IGNORE);
            MPI_Recv(currentCorruptedBlock.hu.getRawPointer(),
                    dataArraySize, MPI_FLOAT, reloadReplica,
                    MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                    interTeamComm, MPI_STATUS_IGNORE);
            std::cout << "T" << myTeam << "R" << myRankInTeam
                    << " : b,h,hv,hu for my block " << i
                    << " are received! Thanks replica " << reloadReplica
                    << std::endl;
            /* compute and validate again */
            currentCorruptedBlock.computeNumericalFluxes();
            bool admissible = currentCorruptedBlock.validateAdmissibility(t);
            if (!admissible) assert(false); // TODO we must abort right ?
            /* problem solved for this corrupted block */
            else {
                primaryBlocksCorrupted[i] = 0;
                std::cout << "T" << myTeam << "R" << myRankInTeam
                        << " : problem solved for my primaryBlock" << i
                        << std::endl;
                //MPI_Send(primaryBlocksCorrupted+i, 1, MPI_BYTE, reloadReplica, 100, interTeamComm); TODO we assume we solved the SDC
            }
        }
    }
}

void tools::Reports::receivePrimaryBlocksReport() {
    *SDC_inReplica = 0;
    /* receive report from other replicas */
    for (int sourceTeam = 0; sourceTeam < numTeams; sourceTeam++) {
        if (sourceTeam != myTeam) {
            MPI_Recv(replicaCorrupted+sourceTeam, 1, MPI_BYTE, sourceTeam,
                     MPI_TAG_REPORT_PRIMARY_BLOCK, interTeamComm, MPI_STATUS_IGNORE);
            if (replicaCorrupted[sourceTeam] == 1) *SDC_inReplica = 1;
        }
    }
}

bool tools::Reports::isLowestHealthyReplica() {
    bool lowestHealthyReplica = true;
    for (unsigned int i = 0; i < myTeam; i++) {
        lowestHealthyReplica &= replicaCorrupted[i];
    }
    return lowestHealthyReplica;
}

void tools::Reports::sendReloadReplica() {
    /* iterate the corrupted replicas */
    for (unsigned int replica = 0; replica < numTeams; replica ++) {
        /* if this replica has reported SDC */
        if (replicaCorrupted[replica] == 1) {
            /* send the reload replica by tag 23 */
            MPI_Send(&myTeam, 1, MPI_INT, replica,
                     MPI_TAG_RECEIVE_RELOAD_REPLICA, interTeamComm);
        }
    }
}

void tools::Reports::recoverCorruptedReplicas() {
    unsigned char secondaryBlocksCorrupted[blocksPerRank - decompFactor];
    /* iterate through all the blocks */
    for (unsigned int i = decompFactor; i < blocksPerRank; i++) {
        // Get a reference to the current block
        const int& currentBlockNr = myBlockOrder[i];
        int currentReplica = currentBlockNr % numTeams;
        /* if this replica has reported SDC */
        if (replicaCorrupted[currentReplica] == 1) {
            /* receive report for the block */
            MPI_Recv(secondaryBlocksCorrupted+(i-decompFactor), 1, MPI_BYTE,
                     currentReplica, MPI_TAG_REPORT_PRIMARY_BLOCK,
                     interTeamComm, MPI_STATUS_IGNORE);
            /* send if it is corrupted */
            if (secondaryBlocksCorrupted[i-decompFactor] == 1) {
                auto& currentSecondaryBlock = *simulationBlocks[currentBlockNr];
                /* Size of the arrays b,h,hv,hu */
                const int dataArraySize = currentSecondaryBlock.dataArraySize;
                int source_rank = currentBlockNr % numTeams;
                std::cout << "T" << myTeam << "R" << myRankInTeam
                          << " : sending to replica in team " << source_rank
                          << " the block " << currentBlockNr
                          << std::endl;
                /* send b,h,hv,hu*/
                MPI_Send(currentSecondaryBlock.b.getRawPointer(),
                         dataArraySize, MPI_FLOAT, source_rank,
                         MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                         interTeamComm);
                MPI_Send(currentSecondaryBlock.h.getRawPointer(),
                         dataArraySize, MPI_FLOAT, source_rank,
                         MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                         interTeamComm);
                MPI_Send(currentSecondaryBlock.hv.getRawPointer(),
                         dataArraySize, MPI_FLOAT, source_rank,
                         MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                         interTeamComm);
                MPI_Send(currentSecondaryBlock.hu.getRawPointer(),
                         dataArraySize, MPI_FLOAT, source_rank,
                         MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                         interTeamComm);
                std::cout << "T" << myTeam << "R" << myRankInTeam
                          << " : sent to replica in team " << source_rank
                          << std::endl;
                secondaryBlocksCorrupted[i-decompFactor] = 0;
            }
        }
    }
}

void tools::Reports::recoverCorruptedReplicas_redundant() {
    /* indicates if a block is corrupted in replica */
    unsigned char replicaBlocksCorrupted[decompFactor];
    /* iterate through all the blocks */
    for (unsigned int i = 0; i < decompFactor; i++) {
        for (int currentReplica = 0; currentReplica < numTeams; currentReplica++) {
            if (currentReplica != myTeam) {
                /* if this replica has reported SDC */
                if (replicaCorrupted[currentReplica] == 1) {
                    /* receive report for the block */
                    MPI_Recv(replicaBlocksCorrupted+i, 1, MPI_BYTE,
                             currentReplica, MPI_TAG_REPORT_PRIMARY_BLOCK,
                             interTeamComm, MPI_STATUS_IGNORE);
                    /* send if it is corrupted */
                    if (replicaBlocksCorrupted[i] == 1) {
                        auto& currentBlock = *simulationBlocks[i];
                        /* Size of the arrays b,h,hv,hu */
                        const int dataArraySize = currentBlock.dataArraySize;
                        std::cout << "T" << myTeam << "R" << myRankInTeam
                                  << " : sending to replica in team " << currentReplica
                                  << std::endl;
                        /* send b,h,hv,hu*/
                        MPI_Send(currentBlock.b.getRawPointer(),
                                 dataArraySize, MPI_FLOAT, currentReplica,
                                 MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                                 interTeamComm);
                        MPI_Send(currentBlock.h.getRawPointer(),
                                 dataArraySize, MPI_FLOAT, currentReplica,
                                 MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                                 interTeamComm);
                        MPI_Send(currentBlock.hv.getRawPointer(),
                                 dataArraySize, MPI_FLOAT, currentReplica,
                                 MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                                 interTeamComm);
                        MPI_Send(currentBlock.hu.getRawPointer(),
                                 dataArraySize, MPI_FLOAT, currentReplica,
                                 MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                                 interTeamComm);
                        std::cout << "T" << myTeam << "R" << myRankInTeam
                                  << " : sent to replica in team " << currentReplica
                                  << std::endl;
                        replicaBlocksCorrupted[i] = 0;
                    }
                }
            }
        }
    }
}

void tools::Reports::reportOwners() {
    /* report to owner replicas of the blocks */
    for (int destTeam = 0; destTeam < numTeams; destTeam++) {
        if (destTeam != myTeam) {
            unsigned char report = 0;
            for (unsigned int i = decompFactor; i < blocksPerRank; i++) {
                /* this block is from my replica in the destTeam */
                if (myBlockOrder[i] % numTeams == destTeam) {
                    if (receivedBlocksCorrupted[i] == 1) report = 1;
                }
            }
            receivedBlockReports[destTeam] = report;
            // receivedBlock validation
            MPI_Send(receivedBlockReports+destTeam, 1, MPI_BYTE, destTeam,
                     MPI_TAG_REPORT_RECEIVED_BLOCK, interTeamComm);
        }
    }
}

void tools::Reports::reportSecondaryBlocks() {
    /* iterate all the replicas that we have received */
    for (int destTeam = 0; destTeam < numTeams; destTeam++) {
        /* if blocks received from replica destTeam is faulty */
        if (receivedBlockReports[destTeam] == 1 && destTeam != numTeams) {
            /* send reports for all received blocks from destTeam */
            for (unsigned int i = decompFactor; i < blocksPerRank; i++) {
                if (myBlockOrder[i] % numTeams == destTeam) {
                    MPI_Send(receivedBlocksCorrupted+i, 1, MPI_BYTE,
                             destTeam, MPI_TAG_REPORT_RECEIVED_BLOCK, interTeamComm);
                }
            }
        }
    }
}

void tools::Reports::recoverCorruptedSecondaryBlocks(float t) {
    /* for each corrupted block receive b,h,hv,hu */
    for (unsigned int i = decompFactor; i < blocksPerRank; i++) {
        if (receivedBlocksCorrupted[i] == 1) {
            /* Get a reference to the current corrupted block */
            const int& currentBlockNr = myBlockOrder[i];
            int reloadReplica = currentBlockNr % numTeams;
            auto& currentCorruptedBlock = *simulationBlocks[currentBlockNr];
            /* Size of the arrays b,h,hv,hu */
            const int dataArraySize = currentCorruptedBlock.dataArraySize;
            std::cout << "T" << myTeam << "R" << myRankInTeam
                      << " : receiving b,h,hv,hu for my secondary block " << currentBlockNr
                      << ", block nr. is " << currentBlockNr
                      << std::endl;
            // receivedBlock recovery by tag 22 : receive order is important! --> b,h,hv,hu
            MPI_Recv(currentCorruptedBlock.b.getRawPointer(),
                     dataArraySize, MPI_FLOAT, reloadReplica,
                     MPI_TAG_RECOVERY_RECEIVED_BLOCK,
                     interTeamComm, MPI_STATUS_IGNORE);
            MPI_Recv(currentCorruptedBlock.h.getRawPointer(),
                     dataArraySize, MPI_FLOAT, reloadReplica,
                     MPI_TAG_RECOVERY_RECEIVED_BLOCK,
                     interTeamComm, MPI_STATUS_IGNORE);
            MPI_Recv(currentCorruptedBlock.hv.getRawPointer(),
                     dataArraySize, MPI_FLOAT, reloadReplica,
                     MPI_TAG_RECOVERY_RECEIVED_BLOCK,
                     interTeamComm, MPI_STATUS_IGNORE);
            MPI_Recv(currentCorruptedBlock.hu.getRawPointer(),
                     dataArraySize, MPI_FLOAT, reloadReplica,
                     MPI_TAG_RECOVERY_RECEIVED_BLOCK,
                     interTeamComm, MPI_STATUS_IGNORE);
            std::cout << "T" << myTeam << "R" << myRankInTeam
                      << " : b,h,hv,hu for my secondary Block" << currentBlockNr
                      << " are received! Thanks replica " << reloadReplica
                      << std::endl;
            /* validate again */
            bool admissible = currentCorruptedBlock.validateAdmissibility_dataArrays(t);
            if (!admissible) assert(false); // TODO we must abort right ?
            /* problem solved for this corrupted block */
            else {
                receivedBlocksCorrupted[i] = 0;
                std::cout << "T" << myTeam << "R" << myRankInTeam
                          << " : problem solved for my secondary Block" << currentBlockNr
                          << ", block nr. is " << currentBlockNr
                          << std::endl;
                //MPI_Send(primaryBlocksCorrupted+i, 1, MPI_BYTE, reloadReplica, 100, interTeamComm); TODO we assume we solved the SDC
            }
        }
    }
}

void tools::Reports::receiveSecondaryBlocksReport() {
    *SDC_inReplica = 0;
    /* receive report from other replicas */
    for (int sourceTeam = 0; sourceTeam < numTeams; sourceTeam++) {
        if (sourceTeam != myTeam) {
            MPI_Recv(replicaCorrupted+sourceTeam, 1, MPI_BYTE, sourceTeam,
                     MPI_TAG_REPORT_RECEIVED_BLOCK, interTeamComm, MPI_STATUS_IGNORE);
            if (replicaCorrupted[sourceTeam] == 1) *SDC_inReplica = 1;
        }
    }
}

void tools::Reports::recoverCorruptedReplicas_TS() {
    /* recovery mode, I will send my primary blocks again for the
     * corrupted blocks in the corrupted replica */
    for (int destTeam = 0; destTeam < numTeams; destTeam++) {
        if (replicaCorrupted[destTeam] == 1 && destTeam != myTeam) {
            for (unsigned int i = 0; i < decompFactor; i++) {
                unsigned char corrupted;
                MPI_Recv(&corrupted, 1, MPI_BYTE, destTeam, MPI_TAG_REPORT_RECEIVED_BLOCK,
                         interTeamComm, MPI_STATUS_IGNORE);
                /* if this block is corrupted in replica destTeam
                 * we send it again */
                if (corrupted == 1) {
                    const int& currentBlockNr = myBlockOrder[i];
                    auto& currentPrimaryBlock = *simulationBlocks[currentBlockNr];
                    /* Size of the arrays b,h,hv,hu */
                    const int dataArraySize = currentPrimaryBlock.dataArraySize;
                    std::cout << "T" << myTeam << "R" << myRankInTeam
                              << " : sending to replica in team " << destTeam
                              << " my block " << currentBlockNr
                              << std::endl;
                    /* send b,h,hv,hu*/
                    MPI_Send(currentPrimaryBlock.b.getRawPointer(),
                             dataArraySize, MPI_FLOAT, destTeam,
                             MPI_TAG_RECOVERY_RECEIVED_BLOCK,
                             interTeamComm);
                    MPI_Send(currentPrimaryBlock.h.getRawPointer(),
                             dataArraySize, MPI_FLOAT, destTeam,
                             MPI_TAG_RECOVERY_RECEIVED_BLOCK,
                             interTeamComm);
                    MPI_Send(currentPrimaryBlock.hv.getRawPointer(),
                             dataArraySize, MPI_FLOAT, destTeam,
                             MPI_TAG_RECOVERY_RECEIVED_BLOCK,
                             interTeamComm);
                    MPI_Send(currentPrimaryBlock.hu.getRawPointer(),
                             dataArraySize, MPI_FLOAT, destTeam,
                             MPI_TAG_RECOVERY_RECEIVED_BLOCK,
                             interTeamComm);
                    std::cout << "T" << myTeam << "R" << myRankInTeam
                              << " : sent to replica in team " << destTeam
                              << std::endl;
                }
            }
        }
    }
}
