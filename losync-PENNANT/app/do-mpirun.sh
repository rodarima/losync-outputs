#!/bin/bash
module use /path/to/install/openmpi/module
module load mpi/openmpi-x86_64-gcc550
export LD_LIBRARY_PATH=/path/to/install/tampi/lib:/opt/gcc-5.5.0/lib64:$LD_LIBRARY_PATH
export NANOS6_CONFIG=$PWD/nanos6.toml
mpirun "$@"

