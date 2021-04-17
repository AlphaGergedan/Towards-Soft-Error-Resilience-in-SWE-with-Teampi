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

TODO

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
- now run `./configure --prefix=/where/to/install --with-ft` to enable ULFM
- finally run `make all install`

Then we should setup teaMPI
- TODO

### Compiling ###

TODO

### Running ###

TODO

