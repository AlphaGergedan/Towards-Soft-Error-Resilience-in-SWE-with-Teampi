# Towards Soft Error Resilience in SWE with TeaMPI #
This repository tries to integrate soft error resilience into SWE using the TeaMPI
library. See the original [SWE](https://github.com/TUM-I5/SWE) and
[teaMPI](https://gitlab.lrz.de/hpcsoftware/teaMPI). Integration over
[hard failure tolerance in teaMPI and SWE](https://github.com/xile273/SWE/tree/simon_task_sharing)

## Soft Error Resilience ##
Increased number of components in high performance computing due to demand on multiple
cores and memories may lead to increased silent data corruptions (SDCs) materialized
as bitflips. In case of an SDC, we are unable to detect the error without observing
the application's outputs because the application seems to be running correctly
but it doesn't instantly report an error when an SDC occurs. In order to provide
resilience for SDCs, we introduce the following methods.

### Method 1 | No Resilience (NoRes) ###
This method doesn't have any built in soft error detection or resilience mechanisms,
so it also doesn't have their overheads. This method can be run to see the fastest
possible runtime and it represents our baseline model. A naive soft error
detection is still possible if the application is run twice and their results are
compared. For correction we need to run it again in order to decide which one of
the initial runs had an SDC by employing a voting mechanism (guess).

SWE domain is divided into `number_of_ranks` blocks and each rank computes its
own block. (See
[swe\_noRes.cpp](https://gitlab.lrz.de/AtamertRahma/towards-soft-error-resilience-in-swe-with-teampi/-/blob/master/src/tolerance/swe_noRes.cpp)
for the implementation and [its executable](#running) for usage examples)

We try to improve this naive method with the following techniques.

### Method 2 | Soft Error Detection Using Hashes (Hashes) ###
In this method we provide soft error detection by employing process replication using 2
TeaMPI teams. We hash the 'tasks' computed by the replicas and share the hashes
across them. This can provide an early SDC detection, however
it cannot correct the soft errors. For soft error resilience an additional run is
required, or this method can be extended to work with 3 Teams again with a voting
mechanism.

We use `std::hash` for better performance. According to
[this stackoverflow page](https://stackoverflow.com/questions/19411742/what-is-the-default-hash-function-used-in-c-stdunordered-map)
the function object `std::hash<>` is based on ***MurmurHashUnaligned2*** for the
specialization for `std::strings`, and it is 8 bytes long. Other hash functions
like SHA-1 could also be used but it is expected to be slower due to its
greater length. We have prepared a Hasher class in `src/tools/Hasher.cpp` for hashing the
data structures of SWE.

We compute one block per rank in team. So we divide the SWE domain into
`number_of_ranks_in_a_team` blocks. Each rank computes its own block and issues
a single heartbeat in each `heartbeat-interval` in simulation time. In each
single heartbeat they send a combined hash value (combined using xor operations)
of the hash values of the updated data arrays and the agreed time step size to their replicas.
TeaMPI handles the comparison of the hashes with our additional modifications, and
makes an early detection of an SDC possible. This method however doesn't guarantee
any hard failure resilience as the heartbeat intervals are bound to simulation time
and not on wall-clock time which is neccessary for detection of any failing or slowing
ranks.
(See
[swe\_softRes\_hashes.cpp](https://gitlab.lrz.de/AtamertRahma/towards-soft-error-resilience-in-swe-with-teampi/-/blob/master/src/tolerance/swe_softRes_hashes.cpp)
for the implementation and [its executable](#running) for usage examples)

We assume that any bitflip that occurs in our data that we hash in each iteration will
eventually be detected by the hash comparisons. This allows us to also detect any SDC
that may have occurred in other data structures if they eventually affect our data
arrays that we are hashing. One downside of this method is that it uses replication
which loses half of our computation power, and even if an SDC is detected we must
restart the application again for recalculation since we don't have a 3 Team support for
this method yet. This can be done again by using more than 2 teams and with a voting
mechanism to identify healthy teams. A healthy team can then write reactive
checkpoints for the corrupted teams and the application can continue its execution.

### Method 3 | Soft Error Resilience Using Admissibility Checks and Task Sharing (Sharing) ###
In order to reduce the overhead of the redundant computation Simon Schuck has integrated
task sharing into SWE for hard error resilience in his work
[here](https://github.com/xile273/SWE/tree/simon_task_sharing). However, task
sharing cannot be integrated if the method depends on redundant computation like
in "Hashes". Therefore we need to use another
soft error detection technique: We validate the results of task computations by
prechecking some predefined admissibility criteria. We use the following 2 main
criteria:
1. ***Physical Admissibility*** : Computations are within some certain physical
constraints, like no negative water height or bathymetry data is constant. These
criteria can be are higly application specific
2. ***Numerical Admissibility*** : No floating point errors are present (NaN)
and relaxed discrete maximum principle (rDMP) in the sense of polynomials.

We divide the global ranks into two teams, and we divide the simulation domain into
multiple blocks, where each rank can have ***primary*** and ***secondary*** blocks.
Each rank will have `d` primary blocks, and
`'number_of_teams - 1' * d` (also equals to number of its replicas) secondary blocks
computed by their replicas, where `d` is the ***decomposition factor***. Ranks will
first compute their primary blocks and check their solutions using some the following
simple predefined admissibility criteria:
1. ***Physical Admissibility Criteria***
    - bathymetry data must be constant
    - water height cannot be negative
2. ***Numerical Admissibility Criteria***
    - no floating point errors (NaN)
    - rDMP

All the criteria except rDMP can be checked cheaply with library functions. For
the implementation see the function ***validateAdmissibility*** in the file
[DimSplitMPIOverdecomp.cpp](https://gitlab.lrz.de/AtamertRahma/towards-soft-error-resilience-in-swe-with-teampi/-/blob/master/src/blocks/DimSplitMPIOverdecomp.cpp).

After validating their primary blocks, ranks share their primary blocks to their
replicas only if the results are admissible. We use the shared tasks
immediately after validating them with the admissibility criteria. If
an admissibility criterion fails, we assume that we have detected an SDC. In
that case we receive the data arrays of the possibly corrupted blocks from a possibly
healthy replica. To check the state of the blocks, we send small reports (MPI\_BYTE) to the
replicas. (See
[swe\_softRes\_admiss\_useShared.cpp](https://gitlab.lrz.de/AtamertRahma/towards-soft-error-resilience-in-swe-with-teampi/-/blob/master/src/tolerance/swe_softRes_admiss_useShared.cpp)
for the implementation and [its executable](#running) for usage examples)

Blocks need to be large enough to make the MPI communications in the task sharing
compensate for their calculations. Task sharing was reported as not flexible
in Schuck's thesis, and we have changed the task sharing so that we only share the
data arrays h, hv and hu, which is a lot easier than sharing the net update arrays.
We again use non-blocking MPI communications in the task sharing. If the communication
is slow due to very large message sizes, increasing the decomposition factor may
help.

With the task sharing we can save redundant computation time but this also comes
with a downside. There is still a chance for an SDC, even if the
results are admissible. This means there is a chance of an undetected SDC, which
spreads immediately with the task sharing because the receiving replicas use the
same criteria to validate the received block and they also would not detect the SDC.
In order to cope with this issue, we introduce the following redundant method.

### Method 4 | Soft Error Resilience Using Admissibility Checks and Redundant Computation (Redundant) ###

This method is very similar to our method ***Sharing***. The only difference is that
we only share a task if it is already corrupted in another replica to recover the
replica. We use the same admissibility criteria and reporting system to validate
the blocks, but since we don't share the tasks and not use them immediately, we can
prevent any possible undetected SDC from spreading and corrupting the other
replicas. We assume that an undetected SDC is eventually going to violate some
of the admissibility criterion like rDMP in the later iterations. With this
method we cannot save any computation overhead caused by the
replication like we did in the previous method, but we can provide additional resilience for
later detected SDCs because each replica keeps its results for itself. (See
[swe\_softRes\_admiss\_redundant.cpp](https://gitlab.lrz.de/AtamertRahma/towards-soft-error-resilience-in-swe-with-teampi/-/blob/master/src/tolerance/swe_softRes_admiss_redundant.cpp)
for the implementation and [its executable](#running) for usage examples)

## Hard Error Resilience ###

We did not test if the hard error resilience is still working in this work.
However the proposed soft error resilient methods are prepared for future hard
resilience integration. For example the methods ***Sharing*** and ***Redundant***
are sending heartbeats using teaMPI in regular intervals that are bound to wall-clock
time. This means we can eventually detect any slowing or failing ranks.

We have also prepared (not tested since hard error resilience was out of scope for this work)
another version of the method ***Hashes*** in the file
[swe\_softRes\_hardRes\_hashes.cpp](https://gitlab.lrz.de/AtamertRahma/towards-soft-error-resilience-in-swe-with-teampi/-/blob/master/src/tolerance/swe_softRes_hardRes_hashes.cpp)
where we send two independent heartbeats between the replicas. One heartbeat
pair is bound to a wall-clock interval that the user specifies for hard error
resilience, and the other single heartbeat is for the soft error
detection which is only bound to the application's simulation time and its only
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
- run `cmake -DMPI_HOME=/ulfm/installation/path -DSOLVER=hlle -B build-directory -S .`
(optional: -DSOLVER={hlle,fwave,augrie} specifies the wave propagation solver,
-DUSE\_DEBUG='on' option can be used for debug mode for developers and -DUSE\_PROFILING='on'
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

Also read the usage outputs of the each executables before executing them by
running `./executable -h`

We have the following executables:
1. ***NoRes***: `./build-directory/swe_noRes` with options
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
2. ***Hashes*** : `./build-directory/swe_softRes_hashes` with options
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
3. ***Sharing*** : `./build-directory/swe_softRes_admiss_useShared`
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
4. ***Redundant*** : `./build-directory/swe_softRes_admiss_redundant`
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

### Testing ###

For testing we have prepared a bash script `runTests.sh` that runs all
the methods and compares their results without injecting any SDC. Make sure
to install [bc](https://www.gnu.org/software/bc/) and export the environment
variable `TEAMS=2` before running this script. One can specify the simulation
parameters like simulation duration or decompositon factor, which only affects the
method ***Sharing*** and ***Redundant***. We take the outputs of the method ***NoRes***
as a reference solution for all the comparisons. We have noticed that netcdf outputs
can sometimes produce different binaries even for the same method for large simulation grids.
In this case one can verify the data arrays' equality by running the script `scripts/ncfloatcmp.py`.
The parameters of the test script are almost
identical to the simulation, the user should also provide a build directory.
Here is an example test run:
```
bash runTests.sh -b build-directory -n 4 -x 200 -y 200 -t 10 -d 8
```

We have also prepared a script to analyze the outcome of SDCs in the methods
***Sharing*** and ***Redundant***. It is important to mention that the
SDC injection is called in the application itself, and therefore can only be changed
in the methods (bitflips are implemented in the blocks with different types
of injections) to test different kinds of SDCs. Default bitflip injection covers the
largest data structures in SWE applications: 12 arrays with
all the main data arrays b,h,hv,hu and the net-update arrays. The array, float
and bit to corrupt is randomly selected. After the injection we wait for the
application to finish and we compare its teams' outputs with a reference
solution from ***NoRes***. We assume that a method has detected and corrected the error, if any
of its team outputs are matching with the reference solution and the application
actually reported the faulty team. So the user can know which teams have detected
an SDC in their blocks and could have written a faulty output. For more detailed
explanation on the outcome analysis, please see the script `scripts/extractSDC_outcomeRate.py`
which was used to analyze the standard output of our bitflip injection script `runSDCAnalysis.sh`
to decide whether an injected error has resulted in a correctable outcome, DUE or SDC. The script
`scripts/runSDCAnalysis.sh` runs the method multiple times and injects a bitflip in
each run. The number of runs with SDC injection can be specified in the `-r` flag.
Provide big numbers to see more precise DUE/Correctable outcome rates. Here is an
example run with 5 SDC injections executions the methods 5 times and injects a random SDC in each of them:
```
bash runSDCAnalysis.sh -b build-directory -n 2 -x 200 -y 200 -t 20 -d 1 -r 5 >> output
python scripts/extractSDC_outcomeRate.py output build-directory 2 200 200 20 1
```
