/**
 * @file src/tools/ftLogger.cpp

 * @author Atamert Rahma
 */

#include "tools/ftLogger.hpp"
#include <string>


/* Public */

tools::FtLogger::FtLogger(unsigned int myTeam, unsigned int myRankInTeam) {
    FtLogger::myTeam = myTeam;
    FtLogger::myRankInTeam = myRankInTeam;

    // get hostname of the machine
    gethostname(FtLogger::hostname, HOST_NAME_MAX);

    // set the strings TODO
    boxInfoTeam = "\u25E6 team: ";
    boxInfoTeam.append(std::to_string(myTeam));
    boxInfoRank = "\u25E6 rank: ";
    boxInfoRank.append(std::to_string(myRankInTeam));
    boxInfoTime = "\u25E6 at t = ";
}

void tools::FtLogger::ft_print_spawnStatus() {
    print_team_rank_long();
    std::cout << " spawned at " << hostname << std::endl;
}

void tools::FtLogger::ft_print_HBstart(double sendWallTime, float t) {
    std::string text = "Heartbeat Start";
    print_Box(text, t);
    std::cout << std::endl;

    // TODO find a way to print out sendWallTime
    /*
    std::cout << "\n\t\u250C" << "\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500"
                << "HEARTBEAT START ++++++++++++++\n\t"
                << "  - team:\t\t" << myTeam << "\n\t"
                << "  - rank:\t\t" << myRankInTeam << "\n\t"
                << "at t = " << t << ", at MPI_Wtime() : "
                << sendWallTime << std::endl;
    */
}

void tools::FtLogger::ft_print_HBend(double timeSinceLastHeartbeat, float t) {
    std::string text = "Heartbeat End";
    print_Box(text, t);
    std::cout << std::endl;
    // TODO find a way to print out sendWallTime

    /*
    std::cout << "\n\t+++++++++ "
                << "HEARTBEAT END ++++++++\n\t"
                << "  - team:\t\t" << myTeam << "\n\t"
                << "  - rank:\t\t" << myRankInTeam << "\n\t"
                << " at t = " << t << ", ended after : "
                << timeSinceLastHeartbeat << " wall time" << std::endl;
    */
}

void tools::FtLogger::ft_calculatingTask(float t) {
    std::cout << "T" << myTeam << "R" << myRankInTeam
                << "\tcalculating.. t = " << t
                << std::endl;
}

void tools::FtLogger::ft_writingTimeStep(float t) {
    std::cout << "T" << myTeam << "R" << myRankInTeam
                << "\t\t--> Writing timestep at: " << t
                << std::endl;
}

void tools::FtLogger::ft_print_loop(double timeSinceLastHeartbeat,
                    clock_t heartbeatInterval, float simulationDuration, float t) {
    std::cout << "T" << myTeam << "R" << myRankInTeam
                << "\t--- timeSinceLastHeartbeat = " << timeSinceLastHeartbeat
                << " < heartbeatInterval = " << heartbeatInterval << " | "
                << "t = " << t << " < simulationDuration = " << simulationDuration
                << std::endl;
}

void tools::FtLogger::ft_block_calculatingTask(unsigned int currentBlockNr, float maxTimestep) {
    std::cout << "T" << myTeam << "R" << myRankInTeam
                << " : Block " << currentBlockNr << ", calculated timestep "
                << maxTimestep << std::endl;
}

void tools::FtLogger::ft_block_sending(float t, unsigned int destTeam) {
    std::cout << "T" << myTeam << "R" << myRankInTeam
                << " : Sending t = " << t << " to Team " << destTeam
                << std::endl;
}
void tools::FtLogger::ft_block_received(float recv_t, unsigned int source_rank) {
    std::cout << "T" << myTeam << "R" << myRankInTeam
                << " : Received t = " << recv_t << " from Team " << source_rank
                << std::endl;
}

void tools::FtLogger::ft_SDC_notDetected() {
    std::cout << "-- T" << myTeam << "R" << myRankInTeam
                << " : No SDC is detected" << std::endl;
}

void tools::FtLogger::ft_SDC_detected() {
    //TODO
}

void tools::FtLogger::ft_SDC_fixed() {
    //TODO
}

void tools::FtLogger::ft_SDC_cannotBeFixed() {
    //TODO
}



/* Private */

void tools::FtLogger::print_team_rank_short() {
    std::cout << "T" << myTeam << "R" << myRankInTeam << sep;
}

void tools::FtLogger::print_team_rank_long() {
    std::cout << "Rank " << myRankInTeam << " of Team " << myTeam << sep;
}

void tools::FtLogger::print_extraLine() {
    std::cout << "\n\n\n";
}

void tools::FtLogger::print_Box(std::string text, float t) {
    std::string currentBoxInfoTime = boxInfoTime;
    currentBoxInfoTime.append(std::to_string(t));

    size_t box_size = std::max(text.length() + 2, boxInfoRank.length() + 3);
    box_size = std::max(box_size, currentBoxInfoTime.length() + 3);

    print_extraLine();

    /* Begin the box */
    std::cout << "\t\t" << box_edge_upperLeft;
    for (size_t i = 0; i < box_size; i++) std::cout << box_line_row;
    std::cout << box_edge_upperRight;

    /* Print text in the box */
    std::cout << "\n\t\t" << box_line_col << " " << text;
    for (size_t i = 0; i < (box_size - (text.length() + 1)); i++) std::cout << " ";
    std::cout << box_line_col;

    /* Print team in the box */
    std::cout << "\n\t\t" << box_line_col << "  " << boxInfoTeam;
    for (size_t i = 0; i < (box_size - boxInfoTeam.length()); i++) std::cout << " ";
    std::cout << box_line_col;

    /* Print rank in the box */
    std::cout << "\n\t\t" << box_line_col << "  " << boxInfoRank;
    for (size_t i = 0; i < (box_size - boxInfoRank.length()); i++) std::cout << " ";
    std::cout << box_line_col;

    /* Print time in the box */
    std::cout << "\n\t\t" << box_line_col << "  " << currentBoxInfoTime;
    for (size_t i = 0; i < (box_size - currentBoxInfoTime.length()); i++) std::cout << " ";
    std::cout << box_line_col;

    /* End the box */
    std::cout << "\n\t\t" << box_edge_downLeft;
    for (size_t i = 0; i < box_size; i++) std::cout << box_line_row;
    std::cout << box_edge_downRight;
}
