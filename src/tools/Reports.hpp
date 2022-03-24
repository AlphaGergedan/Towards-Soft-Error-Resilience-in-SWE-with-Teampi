/**
 * @file src/tools/Reports.hpp
 *
 * @brief Handle the reporting and recovery in the soft error resilience techniques.
 *
 * @author Atamert Rahma rahma@in.tum.de
 */

#ifndef REPORTS_HPP
#define REPORTS_HPP

#include <teaMPI.h>
#include "blocks/DimSplitMPIOverdecomp.hpp"
#include "constants.hpp"

namespace tools {
  class Reports {

    public:

      /**
       * Constructor for Reports. Parameters are needed for other functions especially for the
       * recovery.
       */
      Reports(int numTeams, int myTeam, int myRankInTeam, unsigned int decompFactor, unsigned int blocksPerRank,
              std::vector<std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>> &simulationBlocks,
              std::vector<int> &myBlockOrder, unsigned char *primaryBlocksCorrupted,
              unsigned char *receivedBlocksCorrupted, unsigned char *SDC, unsigned char *SDC_inReplica,
              unsigned char *replicaCorrupted, unsigned char *receivedBlockReports, MPI_Comm interTeamComm);

      Reports(int numTeams, int myTeam, int myRankInTeam, unsigned int decompFactor, unsigned int blocksPerRank,
              std::vector<std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>> &simulationBlocks,
              unsigned char *primaryBlocksCorrupted, unsigned char *SDC,
              unsigned char *SDC_inReplica, unsigned char *replicaCorrupted, MPI_Comm interTeamComm);

      /**
       * Sends SDC reports to all the replicas.
       */
      void reportSDC();

      /**
       * Learns the reload replica for the recovery mode.
       */
      int getReloadReplica();

      /**
       * Sends the corrupted flags of the primary blocks to all the replicas.
       *
       * @param reloadReplica should be identical as the output of getReloadReplica
       */
      void reportPrimaryBlocks(int reloadReplica);

      /**
       * Recover and validate the corrupted blocks.
       *
       * recoverCorruptedPrimaryBlocks_redundant is for the 'redundant' technique.
       *
       * @param reloadReplica should be identical as the output of getReloadReplica
       * @param t should be the current simulation time in the current iteration
       */
      void recoverCorruptedPrimaryBlocks(int reloadReplica, float t);
      void recoverCorruptedPrimaryBlocks_redundant(int reloadReplica, float t); // TODO try to use one function for both of them

      /**
       * Receive SDC report from all the replicas.
       *
       * Sets SDC_inReplica to 1 if any replica has reported a possible SDC.
       */
      void receivePrimaryBlocksReport();

      /**
       * Returns true if this rank is the lowest healthy replica.
       */
      bool isLowestHealthyReplica();


      /**
       * Send reload replica rank to all the corrupted replicas.
       */
      void sendReloadReplica();

      /**
       * Recover the possibly corrupted replicas by sending block information.
       *
       * recoverCorruptedReplicas_redundant is for the 'redundat' technique.
       */
      void recoverCorruptedReplicas();
      void recoverCorruptedReplicas_redundant();

      /**
       * Sends SDC reports to all the secondary block owners.
       */
      void reportOwners();

      /**
       * Sends the corrupted flags of the secondary blocks to all the replicas.
       */
      void reportSecondaryBlocks();

      /**
       * Recover and validate the corrupted secondary blocks.
       *
       * @param t should be the current simulation time in the current iteration
       */
      void recoverCorruptedSecondaryBlocks(float t);

      /**
       * Receive SDC report from all the replicas for secondary blocks.
       *
       * Sets SDC_inReplica to 1 if any replica has reported a possible SDC.
       */
      void receiveSecondaryBlocksReport();

      /**
       * Recover the possibly corrupted replicas by sending block information after
       * the task sharing.
       */
      void recoverCorruptedReplicas_TS();

    private:

      int numTeams, myTeam, myRankInTeam;
      unsigned int decompFactor, blocksPerRank;
      std::vector<std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>> &simulationBlocks;
      std::vector<int> &myBlockOrder;
      unsigned char *primaryBlocksCorrupted, *receivedBlocksCorrupted,
                    *replicaCorrupted, *receivedBlockReports, *SDC, *SDC_inReplica,
                    *secondaryBlocksCorrupted;
      MPI_Comm interTeamComm;

  }; // end of class Reports

} // end of namepace tools


#endif // REPORTS_HPP
