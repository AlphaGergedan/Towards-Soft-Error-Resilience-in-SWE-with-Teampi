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

namespace tools
{
    class FtLogger {

        private:

            /* Team and Rank without line breaks or flush */
            void print_team_rank(unsigned int myTeam, unsigned int myRankInTeam) {
                std::cout << "T" << myTeam << "R" << myRankInTeam;
            }

            /* Print a line break for better readability */
            void print_extraLine() {
                std::cout << "\n\n\n"
                          << std::endl;
            }

        public:

            /** Static fault tolerance printer that all implementations can use */
            static FtLogger ft_logger;

            /**
             * Print status with team, rank and hostname
             *
             * @param myTeam Index of my team
             * @param myRankInTeam Index of my rank in my team
             */
            void ft_print_spawnStatus(unsigned int myTeam, unsigned int myRankInTeam) {
                // TODO move the implementation to private and use again
                char hostname[HOST_NAME_MAX];
                gethostname(hostname, HOST_NAME_MAX);
                std::cout << "Rank " << myRankInTeam << " of Team " << myTeam
                          << " spawned at " << hostname << std::endl;
            }

            /**
             * Print after sending out a heartbeat
             * (with no hashes)
             *
             * @param myTeam Index of my team
             * @param myRankInTeam Index of my rank in my team
             * @param sendWallTime Wall time of the send
             * @param t Current simulation time
             */
            void ft_print_HBstart(unsigned int myTeam, unsigned int myRankInTeam, double sendWallTime, float t) {
                std::cout << "\n\t++++++++++++++ "
                          << "HEARTBEAT START ++++++++++++++\n\t"
                          << "  - team:\t\t" << myTeam << "\n\t"
                          << "  - rank:\t\t" << myRankInTeam << "\n\t"
                          << "at t = " << t << ", at MPI_Wtime() : "
                          << sendWallTime << std::endl;
            }

            /**
             * Print after finishing a heartbeat
             * (with no hashes)
             *
             * @param myTeam Index of my team
             * @param myRankInTeam Index of my rank in my team
             * @param timeSinceLastHeartbeat Wall time between the start and end heartbeat
             * @param t Current simulation time
             */
            void ft_print_HBend(unsigned int myTeam, unsigned int myRankInTeam, double timeSinceLastHeartbeat, float t) {
                std::cout << "\n\t+++++++++ "
                          << "HEARTBEAT END ++++++++\n\t"
                          << "  - team:\t\t" << myTeam << "\n\t"
                          << "  - rank:\t\t" << myRankInTeam << "\n\t"
                          << " at t = " << t << ", ended after : "
                          << timeSinceLastHeartbeat << " wall time" << std::endl;
            }

            /**
             * Call before calculating a timestep
             *
             * @param myTeam Index of my team
             * @param myRankInTeam Index of my rank in my team
             * @param t Current simulation time
             */
            void ft_calculatingTask(unsigned int myTeam, unsigned int myRankInTeam, float t) {
                std::cout << "T" << myTeam << "R" << myRankInTeam
                          << "\tcalculating.. t = " << t
                          << std::endl;
            }

            /**
             * Call before writing a timestep
             *
             * @param myTeam Index of my team
             * @param myRankInTeam Index of my rank in my team
             * @param t Current simulation time
             */
            void ft_writingTimeStep(unsigned int myTeam, unsigned int myRankInTeam, float t) {
                std::cout << "T" << myTeam << "R" << myRankInTeam
                          << "\t\t--> Writing timestep at: " << t
                          << std::endl;
            }

            /**
             * @see tolerance/swe_softRes_and_hardRes_wTaskSharing.cpp
             *
             */
            void ft_print_loop(unsigned int myTeam, unsigned int myRankInTeam, double timeSinceLastHeartbeat,
                               clock_t heartbeatInterval, float simulationDuration, float t) {
                std::cout << "T" << myTeam << "R" << myRankInTeam
                          << "\t--- timeSinceLastHeartbeat = " << timeSinceLastHeartbeat
                          << " < heartbeatInterval = " << heartbeatInterval << " | "
                          << "t = " << t << " < simulationDuration = " << simulationDuration
                          << std::endl;
            }

            /**
             *
             * @see tolerance/swe_softRes_and_hardRes_wTaskSharing.cpp
             */
            void ft_block_calculatingTask(unsigned int myTeam, unsigned int myRankInTeam,
                                          unsigned int currentBlockNr, float maxTimestep) {
                std::cout << "T" << myTeam << "R" << myRankInTeam
                          << " : Block " << currentBlockNr << ", calculated timestep "
                          << maxTimestep << std::endl;
            }

            void ft_block_sending(unsigned int myTeam, unsigned int myRankInTeam,
                                  float t, unsigned int destTeam) {
                std::cout << "T" << myTeam << "R" << myRankInTeam
                          << " : Sending t = " << t << " to Team " << destTeam
                          << std::endl;
            }
            void ft_block_received(unsigned int myTeam, unsigned int myRankInTeam,
                                   float recv_t, unsigned int source_rank) {
                std::cout << "T" << myTeam << "R" << myRankInTeam
                          << " : Received t = " << recv_t << " from Team " << source_rank
                          << std::endl;
            }

            void ft_SDC_notDetected(unsigned int myTeam, unsigned int myRankInTeam) {
                std::cout << "-- T" << myTeam << "R" << myRankInTeam
                          << " : No SDC is detected" << std::endl;
            }

            void ft_SDC_detected() {
                //TODO
            }

            void ft_SDC_fixed() {
                //TODO
            }

            void ft_SDC_cannotBeFixed() {
                //TODO
            }

    }; // end of class FtLogger

} // end of namepace tools

#endif // FT_LOGGER_HPP
