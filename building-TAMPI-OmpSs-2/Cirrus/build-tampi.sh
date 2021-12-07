#!/bin/bash
set -e

# Environment configuration
#module load intel-mpi-18
#module load openmpi/4.1.0
module load mpt
#module load intel-20.4/mpi
module load gcc/8.2.0
module load boost
module load spack
module load autoconf-2.69-gcc-8.2.0-bkc32sr
module load automake-1.16.2-gcc-8.2.0-lyeujjm
module load libtool-2.4.6-gcc-8.2.0-3l4qrtz

export PREFIX=/lustre/home/shared/dc025/ompss2-2021.06-mpt+extrae
export TARGET=${PREFIX}/tampi-install
export TARGET_DEBUG=${PREFIX}/tampi-install-debug
export BOOST_DIR=/lustre/sw/boost/1.73.0
mkdir -p $TARGET
mkdir -p $TARGET_DEBUG
export CC=mpicc CXX=mpicxx FC=mpif90

cd tampi
autoreconf -fiv
#make clean
./configure --prefix=$TARGET --with-boost=$BOOST_DIR
make
make install
cd ..

cd tampi-debug
autoreconf -fiv
#make clean
./configure --prefix=$TARGET_DEBUG --with-boost=$BOOST_DIR --enable-debug-mode
make
make install
cd ..
