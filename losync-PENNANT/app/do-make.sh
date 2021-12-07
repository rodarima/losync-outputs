#!/bin/bash
set -e

INSTALLROOT=/path/to/install

for prog in extrae mercurium nanos6; do
  INSTALLPATH=$INSTALLROOT/$prog
  export PATH=$INSTALLPATH/bin:$PATH
  export LD_LIBRARY_PATH=$INSTALLPATH/lib:$LD_LIBRARY_PATH
  export CFLAGS="-I${INSTALLPATH}/include $CFLAGS"
  export CPPFLAGS=$CFLAGS
  export CXXFLAGS=$CFLAGS
done

export TAMPI_HOME=${INSTALLROOT}/tampi

# Newer GCC fails
export GCC_HOME=/opt/gcc-5.5.0
export PATH=${GCC_HOME}/bin:$PATH

# MPICH
#export PATH=/opt/mpich-gcc-5.5.0/bin:$PATH
#export MPICH_CC=mcc
#export MPICH_CXX=mcxx

# OpenMPI
#module load mpi/openmpi-x86_64
module use ${INSTALLROOT}/openmpi/module
module load mpi/openmpi-x86_64-gcc550
export OMPI_MPICC=mcc
export OMPI_MPICXX=mcxx

make

