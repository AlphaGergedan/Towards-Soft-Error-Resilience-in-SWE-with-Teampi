# Towards Soft Error Resilience in SWE with teaMPI #

This repository tries to integrate soft error resilience in SWE using teaMPI
library. See the original [SWE](https://github.com/TUM-I5/SWE) and
[teaMPI](https://gitlab.lrz.de/hpcsoftware/teaMPI). Integration over
[hard failure tolerance in teaMPI and SWE](https://github.com/xile273/SWE/tree/simon_task_sharing)

## Soft Error Resilience ##

Increased error rates in high performance computing due to multiple cores and
memories may lead to silent data corruptions. In such a case, we are unable to
detect the error without observing the application's results. In order to
provide such resilience, we need to use replication and compare the results of
the replicated processes. A naive method would be to run 2 teams computing the
same application and then compare their results at the end. We try to improve
this by hashing the 'tasks' computed by the replicas and share them across the
replicas. We test 2 different hash functions:
1. The function object std::hash<>, which according to
[this](https://stackoverflow.com/questions/19411742/what-is-the-default-hash-function-used-in-c-stdunordered-map)
page is based on MurmurHashUnaligned2 for the specialization for std::strings.
It is 8 bytes long
2. SHA-1, which is relatively slower but 20 bytes long

We implemented following methods to provide soft error resilience to SWE
applications.

### Single Heartbeats with Hashes ###

This is a soft error resilience technique, where we compute one block per rank
in team. So we divide the SWE domain into `number_of_ranks_in_a_team` blocks.
Each rank computes it's own block and issues a single heartbeat in each
`heartbeat-interval` in simulation time. In each single heartbeat they send the
hashes of the computed updates to their replicas. TeaMPI handles the comparison
of the hashes, and therefore an early detection of a silent data corruption is
possible. This method does not however guarantee any hard failure resilience as
the heartbeat intervals depend on the time steps.
(See
[swe\_softRes.cpp](https://gitlab.lrz.de/AtamertRahma/towards-soft-error-resilience-in-swe-with-teampi/-/blob/master/src/tolerance/swe_softRes.cpp)
for the implementation and [its executable](#running) for usage examples)

### TODO Integrating Single Heratbeats to Hard Failure Resilience ###

One heartbeat pair is bound to the wall clock interval that user specifies for
hard error resilience. And the other heartbeat interval is for the soft error
detection, which only depends on the algorithm and only sending hashes

Theoretically, correction of the data is possible in case of a soft error by
running 3 teams at least, and by deploying a voting mechanism like in
[redMPI](https://www.christian-engelmann.info/?page_id=1873).

## Setting Up ##

First, we need to add User Level Failure Mitigation
[(ulfm2)](https://fault-tolerance.org/2019/11/18/ulfm-4-0-2u1/) to Open MPI
version 4. Follow the instruction in the `INSTALL` file, see the source code
[here](https://bitbucket.org/icldistcomp/ulfm2/src/ulfm/). If the following
steps does not work, read the `HACKING.md`, or try to build using the master
branch of Open MPI itself: https://github.com/open-mpi/ompi
- run `./autogen.pl`
- make sure to have recent versions of GUN Automake and Libtool. Autoconf seems
to have backwards compatibility issues, versions >2.7 may generate faulty
configure file. Version 2.69 seems to be working fine
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
  - `-m HASH_METHOD`: which hash function to use: (0=NONE | 1=stdhash | 2=SHA1),
  default: 1
  - `-c HASH_COUNT`: number of total hashes to send to the replica
  - `-f INJECT_BITFLIP`: injects a bit-flip to the first rank right after the
  simulation time reaches the given argument (double)
  - `-v`: let the simulation produce more output, default: no
  - `-h`: show this help message

  Example run:
  ```
  mpirun --oversubscribe -np 2 build/directory/swe_softRes -t 2 -x 1000 -y 1000 -o outputFile -w -m 1 -c 5 -f 0.2 -v
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

