# Towards Soft Error Resilience in SWE with teaMPI #

This repository tries to integrate soft error resilience in SWE using teaMPI
library. See the original [SWE](https://github.com/TUM-I5/SWE) and
[teaMPI](https://gitlab.lrz.de/hpcsoftware/teaMPI). Integration over
[hard failure tolerance in teaMPI and SWE](https://github.com/xile273/SWE/tree/simon_task_sharing)

## Soft Error Resilience ##

Increased error rates in high performance computing due to multiple cores and
memories may lead to silent data corruptions. In such a case, we are unable to
detect the error without observing the application's results.

Keep in mind:
- don't break teaMPI's asynchronization
- don't break warmSpares. So we don't break fault tolerance

### First Method: Comparing Hashes Using Heartbeats ###

Two possible hashing methods:
  1. SHA-1 hash (20 bytes)
  2. std::hash (size\_t bytes == 8 bytes)

This is a naive soft error resilience technique, where we compute one block per
rank in team. So we divide the domain into 'number of ranks in a team' blocks.
Each rank computes it's own block and issues heartbeats in each `heartbeat-interval`
wall-clock seconds sending the hashes of the computed updates to their replicas.
teaMPI handles the comparison of the hashes, and therefore detection of a silent
data corruption is possible. (See [running examples](#running))

Theoretically, correction of the data is possible in case of a soft error by
running 3 teams at least, and by deploying a voting mechanism like in
[redMPI](https://www.christian-engelmann.info/?page_id=1873).

### Second Method: ###

TODO

## Setting Up ##

First, we need to add User Level Failure Mitigation
[(ulfm2)](https://fault-tolerance.org/2019/11/18/ulfm-4-0-2u1/) to Open MPI
version 4. Follow the instruction in the `INSTALL` file, see the source code
[here](https://bitbucket.org/icldistcomp/ulfm2/src/ulfm/). If the following
steps does not work, read the `HACKING.md`, or try to build using the master
branch of Open MPI itself: https://github.com/open-mpi/ompi
- run `./autogen.pl`
- make sure to have recent versions of GUN Autoconf, Automake, and Libtool
- now run `./configure --prefix=/ulfm/installation/path --with-ft` to enable
ULFM and to install it in the specified `/ulfm/installation/path` directory.
- finally run `make all install` to finish installation.

Also make sure that you have NetCDF library installed.

### Compiling ###

Follow the following steps to compile the project:
- go to the root `/` directory of the project
- make sure you have teaMPI by running: `git submodule update --init --recursive`
- run `cmake -DMPI_HOME=/ulfm/installation/path -B build/directory -S .`
- now prepend the include directory of the installed Open MPI with ULFM (for
example by running `export CPATH=/ulfm/installation/path/include:$CPATH`)
- finally run `make -C build/directory`

Now you should see the executables in the specified 'build/directory'.

### Running ###

You should read the 'README.md' on the submodule before running the simulation.

Also read the 'help' outputs of the each executables before executing them by
running `./executable -h`

An example run for naive soft error detection and hard erpo:
  - `export TEAMS=2` running 2 teams is enough for soft error detection
  - `export SPARES=1` at least one spare process for the hard failure mitigation
  - then run
  ```
  mpirun --oversubscribe -np 2 buid/directory/swe_softRes_and_hardRes_woutTaskSharing -t 10 -x 1000 -y 1000 -i 2 -m 1 -o out/directory/outputName -b out/directory/backupName
  ```

  TODO testing the mitigation methods


```
TODO: add the environment variables to the README on teampi submodule, list:
      SPARES=0                          number of spare processes (warm)
      TEAMS=2                           number of teams
      TMPI_STATS_FILE                   tmpi stats will be written to this file
      TMPI_STATS_OUTPUT_PATH            on this path

      TMPI_FILE
      TMPI_OUTPUT_PATH


 LOGINFO                           enables extra logging information on stdout
                                   but needs to be enabled at compile time and
                                   doesn't compile. Fix this TODO
```

