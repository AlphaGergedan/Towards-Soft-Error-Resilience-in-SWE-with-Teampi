/**
 * @file
 * This file is part of SWE.
 *
 * @author Michael Bader (bader AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Univ.-Prof._Dr._Michael_Bader)
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @section DESCRIPTION
 *
 * Setting of SWE, which uses a wave propagation solver and an artificial or ASAGI scenario on multiple blocks.
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <mpi.h>
#include <string>
#include <vector>
#include <ctime>
#include <csetjmp>

#include "blocks/SWE_Block.hh"

#include "writer/Writer.hh"
#include "reader/Reader.hh"

#ifdef ASAGI
#include "scenarios/SWE_AsagiScenario.hh"
#else
#include "scenarios/SWE_simple_scenarios.hh"
#include "scenarios/SWE_LoadNetCdfScenario.hh"
#endif

#ifdef TEAMPI
#include <teaMPI.h>
#endif

#include "tools/args.hh"
#include "tools/help.hh"
#include "tools/Logger.hh"
#include "tools/ProgressBar.hh"

jmp_buf jumpBuffer;
std::shared_ptr<io::Writer> l_writer;
std::string l_restartBase = "";
std::string l_backupMetadataName, l_baseName, l_backupBase;
std::string l_baseNameInput, l_backupBaseInput;

float l_t;
SWE_Block* l_waveBlock = nullptr;
SWE_Scenario *l_scenario = nullptr;
/**
 * Compute the number of block rows from the total number of processes.
 *
 * The number of rows is determined as the square root of the
 * number of processes, if this is a square number;
 * otherwise, we use the largest number that is smaller than the square
 * root and still a divisor of the number of processes.
 *
 * @param numProcs number of process.
 * @return number of block rows
 */
int computeNumberOfBlockRows(int i_numberOfProcesses) {
  int l_numberOfRows = std::sqrt(i_numberOfProcesses);
  while (i_numberOfProcesses % l_numberOfRows != 0) l_numberOfRows--;
  return l_numberOfRows;
}

// Exchanges the left and right ghost layers.
void exchangeLeftRightGhostLayers( const int i_leftNeighborRank,  SWE_Block1D* o_leftInflow,  SWE_Block1D* i_leftOutflow,
                                   const int i_rightNeighborRank, SWE_Block1D* o_rightInflow, SWE_Block1D* i_rightOutflow,
                                   MPI_Datatype i_mpiCol);

// Exchanges the bottom and top ghist layers.
void exchangeBottomTopGhostLayers( const int i_bottomNeighborRank, SWE_Block1D* o_bottomNeighborInflow, SWE_Block1D* i_bottomNeighborOutflow,
                                   const int i_topNeighborRank,    SWE_Block1D* o_topNeighborInflow,    SWE_Block1D* i_topNeighborOutflow,
                                   const MPI_Datatype i_mpiRow);

#ifdef TEAMPI
void createCheckpointDisk(std::vector<int> failed_teams){
  tools::Logger::logger.printString("Creating Checkpoint");
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  l_writer->writeTimeStep( l_waveBlock->getWaterHeight(),
                            l_waveBlock->getDischarge_hu(),
                            l_waveBlock->getDischarge_hv(),
                            l_t);
  if(rank == 0) l_writer->updateMetadataFile(l_backupMetadataName, l_t, 0);
  l_writer->commitBackup();

  int send = 1;

  for(const auto &i : failed_teams){
    MPI_Send(&send, 1, MPI_INT, TMPI_TeamToWorldRank(rank, i), 0, TMPI_GetWorldComm());
  }
  
}

void loadCheckpointDisk(int reloadTeam){
  tools::Logger::logger.printString("Loading Checkpoint");
  MPI_Status status;
  int recv, rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Recv(&recv, 1, MPI_INT,TMPI_TeamToWorldRank(rank, reloadTeam), 0, TMPI_GetWorldComm(), &status);
  if(!l_waveBlock) delete(l_waveBlock);
  if(!l_scenario) delete(l_scenario);
  l_waveBlock = nullptr;
  l_scenario = nullptr;
  l_restartBase = l_backupBaseInput + "_team_" + std::to_string(reloadTeam);
  std::cout << l_restartBase << std::endl;

  longjmp(jumpBuffer, 1);
}
#endif

/**
 * Main program for the simulation on a single SWE_WavePropagationBlock or SWE_WaveAccumulationBlock.
 */
int main( int argc, char** argv ) {
   //Just to let me attach gdb
    {
    volatile int i = 1;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("PID %d on %s ready for attach\n", getpid(), hostname);
    fflush(stdout);
    while (0 == i)
        sleep(5);
    }
  /**
   * Initialization.
   */
  //! MPI Rank of a process.
  int l_mpiRank;
  //! number of MPI processes.
  int l_numberOfProcesses;

    // check if the necessary command line input parameters are given
  tools::Args args;
  
  args.addOption("grid-size-x", 'x', "Number of cell in x direction");
  args.addOption("grid-size-y", 'y', "Number of cell in y direction");
  args.addOption("output-basepath", 'o', "Output base file name");
  args.addOption("backup-basepath", 'b', "Path of the checkpoint file and metadata used to restart in case of failure");
  args.addOption("restart-basepath", 'r', "Path of the Backup that will be loaded when restarting", tools::Args::Argument::Required, false);
  args.addOption("checkpoint-interval", 'i', "Interval to write checkpoints");
#ifdef ASAGI
  args.addOption("bathymetry-file", 'b', "File containing the bathymetry");
  args.addOption("displacement-file", 'd', "File containing the displacement");
  args.addOption("simul-area-min-x", 0, "Simulation area");
  args.addOption("simul-area-max-x", 0, "Simulation area");
  args.addOption("simul-area-min-y", 0, "Simulation area");
  args.addOption("simul-area-max-y", 0, "Simulation area");
  args.addOption("simul-duration", 0, "Simulation time in seconds");
  #endif

   tools::Args::Result ret = args.parse(argc, argv, true);

  switch (ret)
  {
  case tools::Args::Error:
  	return 1;
  case tools::Args::Help:
	  return 0;
  default:
      break;
  }

  //! total number of grid cell in x- and y-direction.
  int l_nX, l_nY;
  
  clock_t l_checkpointInt;

  // read command line parameters
  l_nX = args.getArgument<int>("grid-size-x");
  l_nY = args.getArgument<int>("grid-size-y");
  l_baseNameInput = args.getArgument<std::string>("output-basepath");
  l_backupBaseInput = args.getArgument<std::string>("backup-basepath");
  l_checkpointInt = args.getArgument<clock_t>("checkpoint-interval");
  if(args.isSet("restart-basepath")){
    l_restartBase = args.getArgument<std::string>("restart-basepath");
  }


  #ifdef TEAMPI
  //! TeaMPI team number
  int l_teamNumber;
  std::function<void(std::vector<int>)> create(createCheckpointDisk);
  std::function<void(int)> load(loadCheckpointDisk);
  TMPI_SetCreateCheckpointCallback(&create);
  TMPI_SetLoadCheckpointCallback(&load);
  TMPI_SetErrorHandlingStrategy(TMPI_WarmSpareErrorHandler);
  #endif

  // initialize MPI
  if(setjmp(jumpBuffer) == 0){
    if ( MPI_Init(&argc,&argv) != MPI_SUCCESS ) {
      std::cerr << "MPI_Init failed." << std::endl;
    }
  } 

  #ifdef TEAMPI
  l_teamNumber = TMPI_GetTeamNumber();
  l_baseName = l_baseNameInput + "_team_" + std::to_string(l_teamNumber);
  l_backupBase = l_backupBaseInput + "_team_" + std::to_string(l_teamNumber);
  #else
  l_backupBase = l_backupBaseInput;
  l_baseName = l_baseNameInput;
  #endif

  // determine local MPI rank
  MPI_Comm_rank(MPI_COMM_WORLD,&l_mpiRank);
  printf("SWE rank %d\n", l_mpiRank);
  // determine total number of processes
  MPI_Comm_size(MPI_COMM_WORLD,&l_numberOfProcesses);

  // initialize a logger for every MPI process
  tools::Logger::logger.setProcessRank(l_mpiRank);

  // print the welcome message
  tools::Logger::logger.printWelcomeMessage();

  // set current wall clock time within the solver
  tools::Logger::logger.initWallClockTime( MPI_Wtime() );
  //print the number of processes
  tools::Logger::logger.printNumberOfProcesses(l_numberOfProcesses);

  //! number of SWE_Blocks in x- and y-direction.
  int l_blocksX, l_blocksY;

  // determine the layout of MPI-ranks: use l_blocksX*l_blocksY grid blocks
  l_blocksY = computeNumberOfBlockRows(l_numberOfProcesses);
  l_blocksX = l_numberOfProcesses/l_blocksY;

  // print information about the grid
  tools::Logger::logger.printNumberOfCells(l_nX, l_nY);
  tools::Logger::logger.printNumberOfBlocks(l_blocksX, l_blocksY);


  //! local position of each MPI process in x- and y-direction.
  int l_blockPositionX, l_blockPositionY;

  // determine local block coordinates of each SWE_Block
  l_blockPositionX = l_mpiRank / l_blocksY;
  l_blockPositionY = l_mpiRank % l_blocksY;

  float l_startTime = 0.0f;


  #ifdef ASAGI
  /*
   * Pixel node registration used [Cartesian grid]
   * Grid file format: nf = GMT netCDF format (float)  (COARDS-compliant)
   * x_min: -500000 x_max: 6500000 x_inc: 500 name: x nx: 14000
   * y_min: -2500000 y_max: 1500000 y_inc: 500 name: y ny: 8000
   * z_min: -6.48760175705 z_max: 16.1780223846 name: z
   * scale_factor: 1 add_offset: 0
   * mean: 0.00217145586762 stdev: 0.245563641735 rms: 0.245573241263
   */

  //simulation area
  float simulationArea[4];
  simulationArea[0] = args.getArgument<float>("simul-area-min-x");
  simulationArea[1] = args.getArgument<float>("simul-area-max-x");
  simulationArea[2] = args.getArgument<float>("simul-area-min-y");
  simulationArea[3] = args.getArgument<float>("simul-area-max-y");

  float simulationDuration = args.getArgument<float>("simul-duration");

  SWE_AsagiScenario l_scenario(args.getArgument<std::string>("bathymetry-file"), args.getArgument<std::string>("displacement-file"),
                               simulationDuration, simulationArea);
  #else
  // create a simple artificial scenario
  if(l_restartBase != ""){
    tools::Logger::logger.printString("Reading");
    io::Reader reader(l_restartBase, l_baseName,l_mpiRank, l_numberOfProcesses, l_blockPositionX, l_blockPositionY);
    l_nX = reader.getGridSizeX();
    l_nY = reader.getGridSizeY();
    l_scenario = reader.getScenario();
    l_startTime = reader.getCurrentTime();
  } else{
    l_scenario = new SWE_SplashingConeScenario();
  }
  #endif

  //! number of grid cells in x- and y-direction per process.
  int l_nXLocal, l_nYLocal;

  //! size of a single cell in x- and y-direction
  float l_dX, l_dY;

  //boundary types and position of the scenario used to init Checkpoint Metadata
  std::vector<BoundaryType> l_boundaryTypes(4);
  std::vector<float> l_boundaryPositions(4);

  for(int i=0; i <= 3; i++){
    l_boundaryTypes[i] = l_scenario->getBoundaryType(static_cast<BoundaryEdge>(i));
    l_boundaryPositions[i] = l_scenario->getBoundaryPos(static_cast<BoundaryEdge>(i));
    
  }

  //end time of the scenario used to init checkpoint metadata
  float l_totalTime = l_scenario->endSimulation();
  // compute local number of cells for each SWE_Block
  l_nXLocal = (l_blockPositionX < l_blocksX-1) ? l_nX/l_blocksX : l_nX - (l_blocksX-1)*(l_nX/l_blocksX);
  l_nYLocal = (l_blockPositionY < l_blocksY-1) ? l_nY/l_blocksY : l_nY - (l_blocksY-1)*(l_nY/l_blocksY);

  // compute the size of a single cell
  l_dX = (l_scenario->getBoundaryPos(BND_RIGHT) - l_scenario->getBoundaryPos(BND_LEFT) )/l_nX;
  l_dY = (l_scenario->getBoundaryPos(BND_TOP) - l_scenario->getBoundaryPos(BND_BOTTOM) )/l_nY;

  // print information about the cell size and local number of cells
  tools::Logger::logger.printCellSize(l_dX, l_dY);
  tools::Logger::logger.printNumberOfCellsPerProcess(l_nXLocal, l_nYLocal);

  //! origin of the simulation domain in x- and y-direction
  float l_originX, l_originY;

  // get the origin from the scenario
  l_originX = l_scenario->getBoundaryPos(BND_LEFT) + l_blockPositionX*l_nXLocal*l_dX;;
  l_originY = l_scenario->getBoundaryPos(BND_BOTTOM) + l_blockPositionY*l_nYLocal*l_dY;

  // create a single wave propagation block
  l_waveBlock = SWE_Block::getBlockInstance(l_nXLocal, l_nYLocal, l_dX, l_dY);
  
  
  // initialize the wave propgation block
  if(l_restartBase != ""){
    l_waveBlock->initScenario(l_originX, l_originY, *l_scenario, true, true);
  } else{
     l_waveBlock->initScenario(l_originX, l_originY, *l_scenario, false, true);
  }
  
  //! time when the simulation ends.
  float l_endSimulation = l_scenario->endSimulation();



  

  /*
   * Connect SWE blocks at boundaries
   */
  // left and right boundaries
  tools::Logger::logger.printString("Connecting SWE blocks at left boundaries.");
  SWE_Block1D* l_leftInflow  = l_waveBlock->grabGhostLayer(BND_LEFT);
  SWE_Block1D* l_leftOutflow = l_waveBlock->registerCopyLayer(BND_LEFT);
  if (l_blockPositionX == 0)
    l_waveBlock->setBoundaryType(BND_LEFT, OUTFLOW);

  tools::Logger::logger.printString("Connecting SWE blocks at right boundaries.");
  SWE_Block1D* l_rightInflow  = l_waveBlock->grabGhostLayer(BND_RIGHT);
  SWE_Block1D* l_rightOutflow = l_waveBlock->registerCopyLayer(BND_RIGHT);
  if (l_blockPositionX == l_blocksX-1)
    l_waveBlock->setBoundaryType(BND_RIGHT, OUTFLOW);

  // bottom and top boundaries
  tools::Logger::logger.printString("Connecting SWE blocks at bottom boundaries.");
  SWE_Block1D* l_bottomInflow  = l_waveBlock->grabGhostLayer(BND_BOTTOM);
  SWE_Block1D* l_bottomOutflow = l_waveBlock->registerCopyLayer(BND_BOTTOM);
  if (l_blockPositionY == 0)
    l_waveBlock->setBoundaryType(BND_BOTTOM, OUTFLOW);

  tools::Logger::logger.printString("Connecting SWE blocks at top boundaries.");
  SWE_Block1D* l_topInflow  = l_waveBlock->grabGhostLayer(BND_TOP);
  SWE_Block1D* l_topOutflow = l_waveBlock->registerCopyLayer(BND_TOP);
  if (l_blockPositionY == l_blocksY-1)
    l_waveBlock->setBoundaryType(BND_TOP, OUTFLOW);

  /*
   * The grid is stored column wise in memory:
   *
   *        ************************** . . . **********
   *        *       *  ny+2 *2(ny+2)*         * (ny+1)*
   *        *  ny+1 * +ny+1 * +ny+1 *         * (ny+2)*
   *        *       *       *       *         * +ny+1 *
   *        ************************** . . . **********
   *        *       *       *       *         *       *
   *        .       .       .       .         .       .
   *        .       .       .       .         .       .
   *        .       .       .       .         .       .
   *        *       *       *       *         *       *
   *        ************************** . . . **********
   *        *       *  ny+2 *2(ny+2)*         * (ny+1)*
   *        *   1   *   +1  *   +1  *         * (ny+2)*
   *        *       *       *       *         *   +1  *
   *        ************************** . . . **********
   *        *       *  ny+2 *2(ny+2)*         * (ny+1)*
   *        *   0   *   +0  *   +0  *         * (ny+2)*
   *        *       *       *       *         *   +0  *
   *        ************************** . . . ***********
   *
   *
   *  -> The stride for a row is ny+2, because we have to jump over a whole column
   *     for every row-element. This holds only in the CPU-version, in CUDA a buffer is implemented.
   *     See SWE_BlockCUDA.hh/.cu for details.
   *  -> The stride for a column is 1, because we can access the elements linear in memory.
   */
  //! MPI row-vector: l_nXLocal+2 blocks, 1 element per block, stride of l_nYLocal+2
  MPI_Datatype l_mpiRow;
  #ifndef CUDA
  MPI_Type_vector(l_nXLocal+2, 1          , l_nYLocal+2, MPI_FLOAT, &l_mpiRow);
  #else
  MPI_Type_vector(1,           l_nXLocal+2, 1          , MPI_FLOAT, &l_mpiRow);
  #endif
  MPI_Type_commit(&l_mpiRow);

  //! MPI row-vector: 1 block, l_nYLocal+2 elements per block, stride of 1
  MPI_Datatype l_mpiCol;
  MPI_Type_vector(1,           l_nYLocal+2, 1,           MPI_FLOAT, &l_mpiCol);
  MPI_Type_commit(&l_mpiCol);

  //! MPI ranks of the neighbors
  int l_leftNeighborRank, l_rightNeighborRank, l_bottomNeighborRank, l_topNeighborRank;

  // compute MPI ranks of the neighbour processes
  l_leftNeighborRank   = (l_blockPositionX > 0) ? l_mpiRank-l_blocksY : MPI_PROC_NULL;
  l_rightNeighborRank  = (l_blockPositionX < l_blocksX-1) ? l_mpiRank+l_blocksY : MPI_PROC_NULL;
  l_bottomNeighborRank = (l_blockPositionY > 0) ? l_mpiRank-1 : MPI_PROC_NULL;
  l_topNeighborRank    = (l_blockPositionY < l_blocksY-1) ? l_mpiRank+1 : MPI_PROC_NULL;

  // print the MPI grid
  tools::Logger::logger.cout() << "neighbors: "
                     << l_leftNeighborRank << " (left), "
                     << l_rightNeighborRank << " (right), "
                     << l_bottomNeighborRank << " (bottom), "
                     << l_topNeighborRank << " (top)" << std::endl;

  // intially exchange ghost and copy layers
  exchangeLeftRightGhostLayers( l_leftNeighborRank,  l_leftInflow,  l_leftOutflow,
                  l_rightNeighborRank, l_rightInflow, l_rightOutflow,
                  l_mpiCol );

  exchangeBottomTopGhostLayers( l_bottomNeighborRank, l_bottomInflow, l_bottomOutflow,
                  l_topNeighborRank,    l_topInflow,    l_topOutflow,
                  l_mpiRow );

  // Init fancy progressbar
  tools::ProgressBar progressBar(l_endSimulation, l_mpiRank);

  // write the output at time zero
  tools::Logger::logger.printOutputTime(0);
  progressBar.update(0.);

  std::string l_fileName = generateBaseFileName(l_baseName,l_blockPositionX,l_blockPositionY);
  std::string l_backupName = generateBaseFileName(l_backupBase, l_blockPositionX, l_blockPositionY);
  l_backupMetadataName = l_backupBase + "_metadata";

  //boundary size of the ghost layers
  io::BoundarySize l_boundarySize = {{1, 1, 1, 1}};
  /*
#ifdef WRITENETCDF
  //construct a NetCdfWriter
  io::NetCdfWriter l_writer( l_fileName,
		  l_waveBlock->getBathymetry(),
		  l_boundarySize,
		  l_nXLocal, l_nYLocal,
		  l_dX, l_dY,
          l_originX, l_originY );
#else
  // Construct a VtkWriter
  io::VtkWriter l_writer( l_fileName,
		  l_waveBlock->getBathymetry(),
		  l_boundarySize,
		  l_nXLocal, l_nYLocal,
		  l_dX, l_dY,
		  l_blockPositionX*l_nXLocal, l_blockPositionY*l_nYLocal );
#endif
  */
  l_writer = io::Writer::createWriterInstance(
          l_fileName,
          l_backupName,
          l_waveBlock->getBathymetry(),
          l_boundarySize,
          l_nXLocal, l_nYLocal,
          l_dX, l_dY,
          l_blockPositionX*l_nXLocal, l_blockPositionY*l_nYLocal,
          l_originX, l_originY,
          0, (args.isSet("restart-basepath"))? true : false); 
      
  tools::Logger::logger.printString(l_backupMetadataName);
  if(l_mpiRank == 0){  
    l_writer->initMetadataFile(l_backupMetadataName,l_totalTime, l_numberOfProcesses, l_nX, l_nY, 0, l_boundaryTypes, l_boundaryPositions);
  }
 

  // Write zero time step
  if(!args.isSet("restart-basepath")){
  l_writer->writeTimeStep( l_waveBlock->getWaterHeight(),
                          l_waveBlock->getDischarge_hu(),
                          l_waveBlock->getDischarge_hv(),
                          l_startTime);
  }
  /**
   * Simulation.
   */
  // print the start message and reset the wall clock time
  progressBar.clear();
  tools::Logger::logger.printStartMessage();
  tools::Logger::logger.initWallClockTime(time(NULL));

  //! simulation time.
  l_t = l_startTime;
  progressBar.update(l_t);

  unsigned int l_iterations = 0;
  
  // loop over checkpoints
  while(l_t <= l_endSimulation) {
    double start = MPI_Wtime();
    double elapsed = 0.0;
    
    // do time steps until next checkpoint is reached
    while( elapsed < l_checkpointInt) {
      //reset CPU-Communication clock
      if(l_t >= l_endSimulation) break;
      tools::Logger::logger.resetClockToCurrentTime("CpuCommunication");
      tools::Logger::logger.resetClockToCurrentTime("Checkpoint");

      // exchange ghost and copy layers
      exchangeLeftRightGhostLayers( l_leftNeighborRank,  l_leftInflow,  l_leftOutflow,
                      l_rightNeighborRank, l_rightInflow, l_rightOutflow,
                      l_mpiCol );

      exchangeBottomTopGhostLayers( l_bottomNeighborRank, l_bottomInflow, l_bottomOutflow,
                      l_topNeighborRank,    l_topInflow,    l_topOutflow,
                      l_mpiRow );

      // reset the cpu clock
      tools::Logger::logger.resetClockToCurrentTime("Cpu");

      // set values in ghost cells
      l_waveBlock->setGhostLayer();

      // compute numerical flux on each edge
      l_waveBlock->computeNumericalFluxes();

      //! maximum allowed time step width within a block.
      float l_maxTimeStepWidth = l_waveBlock->getMaxTimestep();
      
      // update the cpu time in the logger
      tools::Logger::logger.updateTime("Cpu");

      //! maximum allowed time steps of all blocks
      float l_maxTimeStepWidthGlobal;

      // determine smallest time step of all blocks
      MPI_Allreduce(&l_maxTimeStepWidth, &l_maxTimeStepWidthGlobal, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

      // reset the cpu time
      tools::Logger::logger.resetClockToCurrentTime("Cpu");

      // update the cell values
      l_waveBlock->updateUnknowns(l_maxTimeStepWidthGlobal);

      // update the cpu and CPU-communication time in the logger
      
      tools::Logger::logger.updateTime("Cpu");
      tools::Logger::logger.updateTime("CpuCommunication");


      // update simulation time with time step width.
      l_t += l_maxTimeStepWidthGlobal;
      l_iterations++;

      // print the current simulation time
      progressBar.clear();
      tools::Logger::logger.printSimulationTime(l_t);
      progressBar.update(l_t);

      elapsed = MPI_Wtime() - start;
      MPI_Bcast(&elapsed, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      
    }

    // print current simulation time
    progressBar.clear();
    tools::Logger::logger.printOutputTime(l_t);
    progressBar.update(l_t);

    #ifdef TEAMPI
    //check for failures in teaMPI world comm
    MPI_Allreduce(nullptr, nullptr, 0, MPI_INT, MPI_MIN, MPI_COMM_SELF);
    tools::Logger::logger.printString("returned from Heartbeat");
    if(l_t >= l_endSimulation){
      l_writer->writeTimeStep( l_waveBlock->getWaterHeight(),
                            l_waveBlock->getDischarge_hu(),
                            l_waveBlock->getDischarge_hv(),
                            l_t);
      if(l_mpiRank == 0) l_writer->updateMetadataFile(l_backupMetadataName, l_t, 0);
      l_writer->commitBackup();
    }
    #else
    // write output
    double startBackup = MPI_Wtime();
    l_writer->writeTimeStep( l_waveBlock->getWaterHeight(),
                            l_waveBlock->getDischarge_hu(),
                            l_waveBlock->getDischarge_hv(),
                            l_t);
    if(l_mpiRank == 0) l_writer->updateMetadataFile(l_backupMetadataName, l_t, 0);
    l_writer->commitBackup();
    double dur = MPI_Wtime() - startBackup;
    tools::Logger::logger.printString("Checkpoint took: " + std::to_string(dur) + " seconds");
    #endif
  }  

  /**
   * Finalize.
   */
#ifdef ASAGI
  // Free ASAGI resources
  l_scenario.deleteGrids();
#endif

  progressBar.clear();

  // write the statistics message
  tools::Logger::logger.printStatisticsMessage();

  // print the cpu time
  tools::Logger::logger.printTime("Cpu", "CPU time");

  // print CPU + Communication time
  tools::Logger::logger.printTime("CpuCommunication", "CPU + Communication time");

  // print the wall clock time (includes plotting)
  tools::Logger::logger.printWallClockTime(time(NULL));

  // printer iteration counter
  tools::Logger::logger.printIterationsDone(l_iterations);

  // print the finish message
  tools::Logger::logger.printFinishMessage();

  // Dispose of the SWE block!
  delete l_waveBlock;
  delete l_scenario;

  std::cout << l_mpiRank << std::endl;
  // finalize MPI execution
  MPI_Finalize();

  return 0;
}


/**
 * Exchanges the left and right ghost layers with MPI's SendReceive.
 *
 * @param i_leftNeighborRank MPI rank of the  left neighbor.
 * @param o_leftInflow ghost layer, where the left neighbor writes into.
 * @param i_leftOutflow layer where the left neighbor reads from.
 * @param i_rightNeighborRank MPI rank of the right neighbor.
 * @param o_rightInflow ghost layer, where the right neighbor writes into.
 * @param i_rightOutflow layer, where the right neighbor reads form.
 * @param i_mpiCol MPI data type for the vertical gost layers.
 */
void exchangeLeftRightGhostLayers( const int i_leftNeighborRank,  SWE_Block1D* o_leftInflow,  SWE_Block1D* i_leftOutflow,
                                   const int i_rightNeighborRank, SWE_Block1D* o_rightInflow, SWE_Block1D* i_rightOutflow,
                                   MPI_Datatype i_mpiCol) {

  MPI_Status l_status;

  // send to left, receive from the right:
  MPI_Sendrecv( i_leftOutflow->h.elemVector(), 1, i_mpiCol, i_leftNeighborRank,  1,
                o_rightInflow->h.elemVector(), 1, i_mpiCol, i_rightNeighborRank, 1,
                MPI_COMM_WORLD, &l_status );

  MPI_Sendrecv( i_leftOutflow->hu.elemVector(), 1, i_mpiCol, i_leftNeighborRank,  2,
                o_rightInflow->hu.elemVector(), 1, i_mpiCol, i_rightNeighborRank, 2,
                MPI_COMM_WORLD, &l_status );

  MPI_Sendrecv( i_leftOutflow->hv.elemVector(), 1, i_mpiCol, i_leftNeighborRank,  3,
                o_rightInflow->hv.elemVector(), 1, i_mpiCol, i_rightNeighborRank, 3,
                MPI_COMM_WORLD, &l_status );

  // send to right, receive from the left:
  MPI_Sendrecv( i_rightOutflow->h.elemVector(), 1, i_mpiCol, i_rightNeighborRank, 4,
                o_leftInflow->h.elemVector(),   1, i_mpiCol, i_leftNeighborRank,  4,
                MPI_COMM_WORLD, &l_status );

  MPI_Sendrecv( i_rightOutflow->hu.elemVector(), 1, i_mpiCol, i_rightNeighborRank, 5,
                o_leftInflow->hu.elemVector(),   1, i_mpiCol, i_leftNeighborRank,  5,
                MPI_COMM_WORLD, &l_status);

  MPI_Sendrecv( i_rightOutflow->hv.elemVector(), 1, i_mpiCol, i_rightNeighborRank, 6,
                o_leftInflow->hv.elemVector(),   1, i_mpiCol, i_leftNeighborRank,  6,
                MPI_COMM_WORLD, &l_status );

}

/**
 * Exchanges the bottom and top ghost layers with MPI's SendReceive.
 *
 * @param i_bottomNeighborRank MPI rank of the bottom neighbor.
 * @param o_bottomNeighborInflow ghost layer, where the bottom neighbor writes into.
 * @param i_bottomNeighborOutflow host layer, where the bottom neighbor reads from.
 * @param i_topNeighborRank MPI rank of the top neighbor.
 * @param o_topNeighborInflow ghost layer, where the top neighbor writes into.
 * @param i_topNeighborOutflow ghost layer, where the top neighbor reads from.
 * @param i_mpiRow MPI data type for the horizontal ghost layers.
 */
void exchangeBottomTopGhostLayers( const int i_bottomNeighborRank, SWE_Block1D* o_bottomNeighborInflow, SWE_Block1D* i_bottomNeighborOutflow,
                                   const int i_topNeighborRank,    SWE_Block1D* o_topNeighborInflow,    SWE_Block1D* i_topNeighborOutflow,
                                   const MPI_Datatype i_mpiRow) {
  MPI_Status l_status;

  // send to bottom, receive from the top:
  MPI_Sendrecv( i_bottomNeighborOutflow->h.elemVector(), 1, i_mpiRow, i_bottomNeighborRank, 11,
                o_topNeighborInflow->h.elemVector(),     1, i_mpiRow, i_topNeighborRank,11,
                MPI_COMM_WORLD, &l_status );

  MPI_Sendrecv( i_bottomNeighborOutflow->hu.elemVector(), 1, i_mpiRow, i_bottomNeighborRank, 12,
                o_topNeighborInflow->hu.elemVector(),     1, i_mpiRow, i_topNeighborRank,    12,
                MPI_COMM_WORLD, &l_status );

  MPI_Sendrecv( i_bottomNeighborOutflow->hv.elemVector(), 1, i_mpiRow, i_bottomNeighborRank, 13,
                o_topNeighborInflow->hv.elemVector(),     1, i_mpiRow, i_topNeighborRank, 13,
                MPI_COMM_WORLD, &l_status);

  // send to top, receive from the bottom:
  MPI_Sendrecv( i_topNeighborOutflow->h.elemVector(),   1, i_mpiRow, i_topNeighborRank,    14,
                o_bottomNeighborInflow->h.elemVector(), 1, i_mpiRow, i_bottomNeighborRank, 14,
                MPI_COMM_WORLD, &l_status );

  MPI_Sendrecv( i_topNeighborOutflow->hu.elemVector(),   1, i_mpiRow, i_topNeighborRank, 15,
                o_bottomNeighborInflow->hu.elemVector(), 1, i_mpiRow, i_bottomNeighborRank, 15,
                MPI_COMM_WORLD, &l_status );

  MPI_Sendrecv( i_topNeighborOutflow->hv.elemVector(),   1, i_mpiRow, i_topNeighborRank,    16,
                o_bottomNeighborInflow->hv.elemVector(), 1, i_mpiRow, i_bottomNeighborRank, 16,
                MPI_COMM_WORLD, &l_status );

}
