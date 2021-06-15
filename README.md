# Towards Soft Error Resilience in SWE with teaMPI #
This repository tries to integrate soft error resilience in SWE using teaMPI
library. See the original [SWE](https://github.com/TUM-I5/SWE) and
[teaMPI](https://gitlab.lrz.de/hpcsoftware/teaMPI). Integration over
[hard failure tolerance in teaMPI and SWE](https://github.com/xile273/SWE/tree/simon_task_sharing)

## Soft Error Resilience ##
Increased error rates in high performance computing due to multiple cores and
memories may lead to silent data corruptions (SDCs) with bitflips. In such a
case, we are unable to detect the error without observing the application's
results because in case of an SDC the application seems to be running correctly
and it doesn't instantly report an error when an SDC occurs. In order to provide
detection and resilience for SDCs, we introduce the following methods.

### Method 1 | No Resilience ###
This method doesn't have any built in soft detection or resilience mechanism, so
it also doesn't have their overheads. This method can be run to see the fastest
possible computation and it represents our baseline model. A naive soft error
detection can be done by running the application twice. For resilience we need
to run it again in order to decide which one of the initial runs had an SDC by
using voting.

SWE domain is divided into `number_of_ranks` blocks and each rank computes its
own block. (See
[swe\_noRes.cpp](https://gitlab.lrz.de/AtamertRahma/towards-soft-error-resilience-in-swe-with-teampi/-/blob/master/src/tolerance/swe_noRes.cpp)
for the implementation and [its executable](#running) for usage examples)

We try to improve this naive method with the following techniques.

### Method 2 | Soft Error Detection Using Hashes ###
In this method we provide soft error detection by running the application
redundantly twice. We hash the 'tasks' computed by the replicas of 2 Teams and
share the hashes across them. This can provide an early SDC detection, however
it cannot correct the soft errors. For soft resilience an additional run is
required, or this method can be extended to work with 3 Teams with a voting
mechanism.

We use std::hash for better performance. According to
[this stackoverflow page](https://stackoverflow.com/questions/19411742/what-is-the-default-hash-function-used-in-c-stdunordered-map)
the function object std::hash<> is based on MurmurHashUnaligned2 for the
specialization for std::strings, and it is 8 bytes long. Other hash functions
like SHA-1 could also be used but it is expected to be slower due to its
increased length.

We compute one block per rank in team. So we divide the SWE domain into
`number_of_ranks_in_a_team` blocks. Each rank computes its own block and issues
a single heartbeat in each `heartbeat-interval` in simulation time. In each
single heartbeat they send the hashes of the computed updates to their replicas.
TeaMPI handles the comparison of the hashes, and therefore an early detection of
an SDC is possible. This method however does not guarantee any hard failure
resilience as the heartbeat intervals depend on the time steps.
(See
[swe\_softRes\_hashes.cpp](https://gitlab.lrz.de/AtamertRahma/towards-soft-error-resilience-in-swe-with-teampi/-/blob/master/src/tolerance/swe_softRes_hashes.cpp)
for the implementation and [its executable](#running) for usage examples)

We assume that any bitflip occurs in our data will eventually be detected by the
hash comparisons. This allows us to detect any SDC that may occur in our data
arrays that we are hashing. But this method uses replication, which loses half
of our computation power, and even if an SDC is detected we must restart the
application again for recalculation since we don't have a 3 Team support for
this method yet. This can be done in further research with a voting mechanism
to choose the healthy team. It can then write reactive checkpoint for the other
teams and the application can continue.

### Method 3 | Soft Error Resilience Using Admissibility Checks and Task Sharing ###
TODO

soft error resilience with admissibility checks and using shared tasks
Checks for admissibility of the computations (also see validateAdmissibility
in src/blocks/DimSplitMPIOverdecomp.cpp) and only share the results if they
are admissible. If they are not, a healthy replica sends its block data to
the failed replicas. We use shared tasks immediately, and check
them for admissibility in case of SDC during transmission (undetected SDC
spreads immediately, but saves computation time for secondary blocks).

### Method 4 | Soft Error Resilience Using Admissibility Checks and Redundant Computation ###
TODO


### Soft and Hard Error Resilience with Task Sharing ###
or we need to check if our results are
"sane" enough to proceed (like checking for absence of NaNs etc.).
techniques, and still try to provide soft resilient HPC applications using
teaMPI.


We already had hard error resilience using warm spares with task sharing in
[Simon's version](https://github.com/xile273/SWE/tree/simon_task_sharing) of SWE
with teaMPI. Integrating soft error resilience to task sharing using hashes is
not possible in theory, hence hashes need to be validated after redundant
computations. Fortunately we can still validate the task computations by
prechecking admissibility criteria of the task computations, which are usually
not transparent to applications. In this version of SWE, we share a computed
task, only if it is validated using the following 2 criteria:
1. ***Physical Admissibility*** : Computations are within some certain physical
constraints, like no negative water height or bathymetry data is constant. These
criteria are higly application specific
2. ***Numerical Admissibility*** : No floating point errors are present (NaN)
and relaxed discrete maximum principle TODO

If any data array, that can suffer from a bit flip, doesn't meet a predefined
admissibility criterion, then we assume that an SDC is present. Unfortunately
we cannot assume that any SDC occurred in the data will be detected in the
future because of the computations depend on the data, if we want to use task
sharing. This means, we have to check all the data array and share our computed
tasks. If we detect an SDC, the other team creates a checkpoint and the
corrupted team loads it. This however requires some additional synchronization
between and within the teams. See [example runs](#running) for usage.



## Hard Error Resilience ###

Method 2 was integrated into the heartbeats that teaMPI uses in
order to detect hard error

### TODO Integrating Single Heratbeats to Hard Failure Resilience ###

One heartbeat pair is bound to the wall clock interval that user specifies for
hard error resilience. And the other heartbeat interval is for the soft error
detection, which only depends on the algorithm and only sending hashes

Theoretically, correction of the data is possible in case of a soft error by
running 3 teams at least, and by deploying a voting mechanism like in
[redMPI](https://www.christian-engelmann.info/?page_id=1873).


### TODO Integrating Hard Error Resilience Into Other Methods ###


## Setting Up ##

First, we need to add User Level Failure Mitigation
[(ulfm2)](https://fault-tolerance.org/2019/11/18/ulfm-4-0-2u1/) to Open MPI
version 4. Follow the instruction in the `INSTALL` file, see the source code
[here](https://bitbucket.org/icldistcomp/ulfm2/src/ulfm/). If the following
steps does not work, read the `HACKING.md`, or try to build using the master
branch of Open MPI itself: https://github.com/open-mpi/ompi
- run `./autogen.pl`
- make sure to have recent versions of GUN Automake and Libtool. Autoconf seems
to have backwards compatibility issues, versions >2.7 may generate a faulty
configuration file. Version 2.69 seems to be working fine
- now run `./configure --prefix=/ulfm/installation/path --with-ft` to enable
ULFM and to install it in the specified `/ulfm/installation/path` directory.
- finally run `make all install` to finish installation.

Also make sure that you have NetCDF library installed. You may need to set the
environment variable `HDF5_USE_FILE_LOCKING=FALSE` for reading netcdf files
depending on your configurations

### Compiling ###

Follow the following steps to compile the project:
- go to the root `/` directory of the project
- make sure you have teaMPI by running: `git submodule update --init --recursive`
- run `cmake -DMPI_HOME=/ulfm/installation/path -B build/directory -S .`
(-DUSE\_DEBUG='on' option can be used for debug mode)
- now prepend the include directory of the installed Open MPI with ULFM (for
example by running `export CPATH=/ulfm/installation/path/include:$CPATH`)
- finally run `make -C build/directory`

Now you should see the executables in the specified 'build/directory'.

### Running ###

You should read the
'[README.md](https://gitlab.lrz.de/AtamertRahma/teampi-soft-error-resilience/-/blob/master/README.md)'
on the submodule before running the simulation and set the required environment
variables.

Also read the 'help' outputs of the each executables before executing them by
running `./executable -h`

For soft error resilience we have the following executables:
1. single heartbeats with hashes: `./build/directory/swe_softRes` with options
  - `-t SIMULATION_DURATION`: time in seconds to simulate
  - `-x RESOLUTION_X`: number of simulated cells in x-direction
  - `-y RESOLUTION_Y`: number of simulated cells in y-direction
  - `-o OUTPUT_BASEPATH`: output base file name
  - `-w`: write output using netcdf writer to the specified output base file
  - `-m HASH_METHOD`: which hash function to use: (0=NONE | 1=stdhash),
  default: 1
  - `-c HASH_COUNT`: number of total hashes to send to the replica
  - `-f INJECT_BITFLIP`: injects a bit-flip to the first rank right after the
  simulation time reaches the given argument (double)
  - `-v`: let the simulation produce more output, default: no
  - `-h`: show this help message

  Example run:
  ```
  mpirun -np 2 build/directory/swe_softRes -t 2 -x 1000 -y 1000 -o outputFile -w -m 1 -c 5 -f 0.2 -v
  ```
2.  soft and hard error resilience with task sharing:
`./build/directory/swe_softRes_and_hardRes_wTaskSharing` with options
  - `-t SIMULATION_DURATION`: time in seconds to simulate
  - `-i HEARTBEAT_INTERVAL` : wall-clock time in seconds to wait between heartbeats
  - `-x RESOLUTION_X`: number of simulated cells in x-direction
  - `-y RESOLUTION_Y`: number of simulated cells in y-direction
  - `-d DECOMP_FACTOR` : Split each rank into "TEAMS" * "DECOMP\_FACTOR" blocks
  - `-o OUTPUT_BASEPATH`: output base file name
  - `-b BACKUP_BASEPATH`: backup base file name
  - `-r RESTART_BASEPATH`: restart base file name
  - `-w`: write output using netcdf writer to the specified output base file
  - `-f INJECT_BITFLIP`: injects a bit-flip to the first rank right after the
  simulation time reaches the given argument (double)
  - `-k KILL_RANK` : kills the rank 0 of team 0 at the specified simulation time
  - `-v`: let the simulation produce more output, default: no
  - `-h`: show this help message

  Example run:
  ```
  mpirun -np 2 build/directory/swe_softRes_and_hardRes_wTaskSharing -t 2 -i 5 -x 1000 -y 1000 -d 1 -o outputFile -b backupFile -w -v
  ```

```
TODO: add the environment variables to the README on teampi submodule, list:
      SPARES=0                          number of spare processes (warm)
      TEAMS=2                           number of teams
      TMPI_STATS_FILE                   tmpi stats will be written to this file
      TMPI_STATS_OUTPUT_PATH            on this path

      TMPI_FILE
      TMPI_OUTPUT_PATH
```

