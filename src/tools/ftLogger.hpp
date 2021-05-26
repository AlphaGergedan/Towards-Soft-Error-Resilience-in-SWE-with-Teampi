/**
 * @file src/tools/ftLogger.hpp
 * @brief Printer for fault tolerance techniques on stdout
 *
 * TODO differentiate the different ft techniques CR + SOFT detection etc.
 * TODO move public implementations to private
 * TODO figure out how to save team and rank
 *
 * @author Atamert Rahma
 */

#ifndef FT_LOGGER_HPP
#define FT_LOGGER_HPP

#include <iomanip>
#include <iostream>
#include <ctime>
#include <ostream>
#include <unistd.h>
#include <limits.h>

namespace tools {
  class FtLogger {

    private:

        /* TODO testing out this, or maybe just employ progress bar */
        std::string sep = " | ";
        std::string box_edge_upperLeft = "\u256D";
        std::string box_edge_upperRight = "\u256E";
        std::string box_edge_downLeft = "\u2570";
        std::string box_edge_downRight = "\u256F";
        std::string box_line_row = "\u2500";
        std::string box_line_col = "\u2502";

        std::string boxInfoTeam;
        std::string boxInfoRank;
        std::string boxInfoTime;

        /* Some teaMPI Information */
        unsigned int myTeam;
        unsigned int myRankInTeam;
        char hostname[HOST_NAME_MAX];

        /*Short Print Team and Rank */
        void print_team_rank_short();

        /* Long Print Team and Rank */
        void print_team_rank_long();

        /* Print 3 line breaks for better readability */
        void print_extraLine();

        /* Prints small box with text, text should not be very large TODO */
        void print_Box(std::string text, float t);

    public:

        /**
         * TODO Constructor
         */
        FtLogger(unsigned int myTeam, unsigned int myRankInTeam);

        /**
         * Print status with team, rank and hostname
         *
         * @param myTeam Index of my team
         * @param myRankInTeam Index of my rank in my team
         */
        void ft_print_spawnStatus();

        /**
         * Print after sending out a heartbeat
         * (with no hashes)
         *
         * @param myTeam Index of my team
         * @param myRankInTeam Index of my rank in my team
         * @param sendWallTime Wall time of the send
         * @param t Current simulation time
         */
        void ft_print_HBstart(double sendWallTime, float t);

        /**
         * Print after finishing a heartbeat
         * (with no hashes)
         *
         * @param myTeam Index of my team
         * @param myRankInTeam Index of my rank in my team
         * @param timeSinceLastHeartbeat Wall time between the start and end heartbeat
         * @param t Current simulation time
         */
        void ft_print_HBend(double timeSinceLastHeartbeat, float t);

        /**
         * Call before calculating a timestep
         *
         * @param myTeam Index of my team
         * @param myRankInTeam Index of my rank in my team
         * @param t Current simulation time
         */
        void ft_calculatingTask(float t);

        /**
         * Call before writing a timestep
         *
         * @param myTeam Index of my team
         * @param myRankInTeam Index of my rank in my team
         * @param t Current simulation time
         */
        void ft_writingTimeStep(float t);

        /**
         * @see tolerance/swe_softRes_and_hardRes_wTaskSharing.cpp
         *
         */
        void ft_print_loop(double timeSinceLastHeartbeat, clock_t heartbeatInterval,
                           float simulationDuration, float t);

        /**
         *
         * @see tolerance/swe_softRes_and_hardRes_wTaskSharing.cpp
         */
        void ft_block_calculatingTask(unsigned int currentBlockNr, float maxTimestep);

        void ft_block_sending(float t, unsigned int destTeam);


        void ft_block_received(float recv_t, unsigned int source_rank);

        void ft_SDC_notDetected();

        // TODO
        void ft_SDC_detected();

        //TODO
        void ft_SDC_fixed();

        //TODO
        void ft_SDC_cannotBeFixed();

  }; // end of class FtLogger
} // end of namepace tools

#endif // FT_LOGGER_HPP
