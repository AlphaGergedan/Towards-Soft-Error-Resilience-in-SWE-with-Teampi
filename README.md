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

We use `std::hash` for better performance. According to
[this stackoverflow page](https://stackoverflow.com/questions/19411742/what-is-the-default-hash-function-used-in-c-stdunordered-map)
the function object `std::hash<>` is based on ***MurmurHashUnaligned2*** for the
specialization for `std::strings`, and it is 8 bytes long. Other hash functions
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
hash comparisons. This allows us to detect any SDC that may have occurred in our data
arrays that we are hashing. But this method uses replication, which loses half
of our computation power, and even if an SDC is detected we must restart the
application again for recalculation since we don't have a 3 Team support for
this method yet. This can be done in further research with a voting mechanism
to choose the healthy team. It can then write reactive checkpoint for the other
teams and the application can continue.

### Method 3 | Soft Error Resilience Using Admissibility Checks and Task Sharing ###
In order to reduce the overhead of the redundant computation Simon Schuck integrated
task sharing into hard error resilience in his work
[here](https://github.com/xile273/SWE/tree/simon_task_sharing). However task
sharing cannot be integrated if the method depends on redundant computation like
in the previous soft resilience methods. Therefore we need to use another
resilience technique. We can validate the results of task computations by
prechecking some predefined admissibility criteria. We use the following 2 main
criteria:
1. ***Physical Admissibility*** : Computations are within some certain physical
constraints, like no negative water height or bathymetry data is constant. These
criteria are higly application specific
2. ***Numerical Admissibility*** : No floating point errors are present (NaN)
and relaxed discrete maximum principle (DMP) in the sense of polynomials.

We divide the ranks into two teams, and we divide the domain into
multiple blocks, where each rank can have ***primary*** and ***secondary*** blocks.
Each rank will have `d` primary blocks, and
`number_of_teams - 1' * d` secondary blocks computed by their replicas,
where `d` is the ***decomposition factor***. Ranks will first compute their
primary blocks and the solutions are checked for certain admissibility criteria.
We use the following simple admissibility criteria:
1. ***Physical Admissibility Criteria***
    - bathymetry data must be constant
    - water height cannot be negative
2. ***Numerical Admissibility Criteria***
    - no floating point errors (NaN)
    - DMP

All the criteria except DMP can be checked cheaply with library functions. For
the implementation see the function ***validateAdmissibility*** in the file
[DimSplitMPIOverdecomp.cpp](https://gitlab.lrz.de/AtamertRahma/towards-soft-error-resilience-in-swe-with-teampi/-/blob/master/src/blocks/DimSplitMPIOverdecomp.cpp).

After validating their primary blocks, ranks share their primary blocks to their
replicas, only if the results are admissible. We use the shared tasks
immediately after validating them with the admissibility criteria. If
an admissibility criterion fails, we assume that we have detected an SDC. In
that case we receive the data arrays of the corrupted blocks from a healthy
replica. To check the state of the blocks, we send small reports (MPI\_BYTE) to the
replicas. (See
[swe\_softRes\_admiss\_useShared.cpp](https://gitlab.lrz.de/AtamertRahma/towards-soft-error-resilience-in-swe-with-teampi/-/blob/master/src/tolerance/swe_softRes_admiss_useShared.cpp)
for the implementation and [its executable](#running) for usage examples)

Blocks need to be large enough to make task sharing less expensive than
their calculation. task sharing was reported as not flexible in Simon's thesis,
we changed the task sharing so that we only share the data arrays b, h, hv and hu,
which is a lot easier than sharing the net update arrays. We also managed to
improve task sharing so that we check if one block is received and update the
results while waiting for other blocks. We use non-blocking MPI communication
to manage this. If the communication is slow, increasing the decomposition
factor may decrease the runtime.

With the task sharing we can save redundant computation time, but there is also
a downside of this method. There is still a chance for an SDC, even if the
results are admissible. This means there is a chance of an undetected SDC, which
spreads immediately with the task sharing, because the receiving replicas use the
same criteria to validate the received block and they also cannot detect the SDC.
In order cope with this issue, we introduce a similar method.

### Method 4 | Soft Error Resilience Using Admissibility Checks and Redundant Computation ###

This method is very similar to the method 3 above. The only difference is that
we only share a task if it is already corrupted in another replica to save the
replica. We use the same admissibility criteria and reporting system to validate
the blocks, but since we don't share the tasks and not use them immediately, we can
prevent any possible undetected SDC from spreading and corrupting the other
replicas. We assume that an undetected SDC is eventually going to violate some
of the admissibility criterion like DMP in the further timesteps. With this
method we cannot save any computation overhead caused by the
replication like we did in the previous method, but we can provide resilience for
later detected SDCs because each replica keeps its results for itself. (See
[swe\_softRes\_admiss\_redundant.cpp](https://gitlab.lrz.de/AtamertRahma/towards-soft-error-resilience-in-swe-with-teampi/-/blob/master/src/tolerance/swe_softRes_admiss_redundant.cpp)
for the implementation and [its executable](#running) for usage examples)

## Hard Error Resilience ###

We did not test if the hard error resilience is still working in this work.
However all the methods (except No Resilience) are prepared for future hard
resilience integration. For example the method 3 and 4 are sending heartbeats
using teaMPI in regular intervals, which does not depend on the simulation time,
but depend on the wall-clock time. This means we can detect the slowing ranks.
The most challenging part of integrating hard error resilience is the handling
of failed MPI calls in case of a hard failure. ULFM provides us error codes, which
can be checked after every MPI call. However we have a lot of MPI calls,
especially in the method 3, where we share the computed tasks in every timestep.

We also prepared another version of the method 2 (Soft Error Detection Using
Hashes) in the file
[swe\_softRes\_hardRes\_hashes.cpp](https://gitlab.lrz.de/AtamertRahma/towards-soft-error-resilience-in-swe-with-teampi/-/blob/master/src/tolerance/swe_softRes_hardRes_hashes.cpp)
where we send two independent heartbeats between the replicas. One heartbeat
pair is bound to a wall-clock interval that the user specifies for hard error
resilience. And the other single heartbeat is for the soft error
detection, which only depends on the application's simulation time and its only
purpose is to send the hashes of the computed data arrays to the other replicas
for validation. This method can be improved in the future
for hard error resilience integration.

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
- run `cmake -DMPI_HOME=/ulfm/installation/path -B build-directory -S .`
(-DUSE\_DEBUG='on' option can be used for debug mode and -DUSE\_PROFILING='on'
option can be used for profiling with a profiler like VTune)
- now prepend the include directory of the installed Open MPI with ULFM (for
example by running `export CPATH=/ulfm/installation/path/include:$CPATH`)
- finally run `make -C build-directory`

Now you should see the executables in the specified 'build-directory'.

### Running ###

You should read the
'[README.md](https://gitlab.lrz.de/AtamertRahma/teampi-soft-error-resilience/-/blob/master/README.md)'
on the submodule before running the simulation and set the required environment
variables (for this repository important one is the `TEAMS` environment
variable, it is 2 by default).

Also read the 'help' outputs of the each executables before executing them by
running `./executable -h`

For soft error resilience we have the following executables:
1. Method 1 - No Resilience: `./build-directory/swe_noRes` with options
  - `-t SIMULATION_DURATION`: time in seconds to simulate
  - `-x RESOLUTION_X`: number of simulated cells in x-direction
  - `-y RESOLUTION_Y`: number of simulated cells in y-direction
  - `-o OUTPUT_BASEPATH`: output base file name
  - `-w`: write output using netcdf writer to the specified output base file
  - `-v`: let the simulation produce more output, default: no
  - `-h`: show this help message
  Example run (hint: you may need to set the `--oversubscribe` flag if you want
  to start more processes than you have):
  ```
  mpirun -np 2 build-directory/swe_noRes -t 2 -x 1000 -y 1000 -o outputFile -w -v
  ```
2. Method 2 - Soft Error Detection Using Hashes: `./build-directory/swe_softRes_hashes` with options
  - `-t SIMULATION_DURATION`: time in seconds to simulate
  - `-x RESOLUTION_X`: number of simulated cells in x-direction
  - `-y RESOLUTION_Y`: number of simulated cells in y-direction
  - `-o OUTPUT_BASEPATH`: output base file name
  - `-w`: write output using netcdf writer to the specified output base file
  - `-m HASH_METHOD`: which hash function to use: (0=NONE | 1=stdhash),
  default: 1
  - `-c HASH_COUNT`: number of total hashes to send to the replica
  - `-f INJECT_BITFLIP`: Injects a random bit-flip into a random data array in a random team and rank right after the simulation time reaches the given time
  - `-v`: let the simulation produce more output, default: no
  - `-h`: show this help message
  Example run:
  ```
  mpirun -np 4 build-directory/swe_softRes_hashes -t 2 -x 1000 -y 1000 -o outputFile -w -m 1 -c 5 -v
  ```
3. Method 3 - Soft Error Resilience Using Admissibility Checks and Task Sharing : `./build-directory/swe_softRes_admiss_useShared`
  - `-t SIMULATION_DURATION`: time in seconds to simulate
  - `-x RESOLUTION_X`: number of simulated cells in x-direction
  - `-y RESOLUTION_Y`: number of simulated cells in y-direction
  - `-o OUTPUT_BASEPATH`: output base file name
  - `-w`: write output using netcdf writer to the specified output base file
  - `-i HEARTBEAT_INTERVAL`: wall-clock time in seconds to wait between heartbeats
  - `-d DECOMP_FACTOR` : split each rank into `TEAMS * decomp-factor` blocks (for better runtimes in some cases)
  - `-f INJECT_BITFLIP`: Injects a random bit-flip into a random data array in a random team and rank right after the simulation time reaches the given time
  - `-v`: let the simulation produce more output, default: no
  - `-h`: show this help message
  Example run:
  ```
  mpirun -np 2 build-directory/swe_softRes_admiss_useShared -t 2 -x 1000 -y 1000 -o outputFile -w -i 5 -d 1 -v
  ```
4. Method 4 - Soft Error Resilience Using Admissibility Checks and Redundant Computation: `./build-directory/swe_softRes_admiss_redundant`
  - `-t SIMULATION_DURATION`: time in seconds to simulate
  - `-x RESOLUTION_X`: number of simulated cells in x-direction
  - `-y RESOLUTION_Y`: number of simulated cells in y-direction
  - `-o OUTPUT_BASEPATH`: output base file name
  - `-w`: write output using netcdf writer to the specified output base file
  - `-i HEARTBEAT_INTERVAL`: wall-clock time in seconds to wait between heartbeats
  - `-d DECOMP_FACTOR` : split each rank into `decomp-factor` blocks
  - `-f INJECT_BITFLIP`: Injects a random bit-flip into a random data array in a random team and rank right after the simulation time reaches the given time
  - `-v`: let the simulation produce more output, default: no
  - `-h`: show this help message
  Example run:
  ```
  mpirun -np 4 build-directory/swe_softRes_admiss_redundant -t 2 -x 1000 -y 1000 -o outputFile -w -i 5 -d 1 -v
  ```
