#!/bin/bash
set -e

#module load intel-mpi-18
module load mpt
#module load openmpi/4.1.0-ucx-gcc8

module load gcc/8.2.0
export TARGET_EXTRAE=/lustre/home/shared/dc025/ompss2-2021.06-mpt+extrae/extrae

# intel-mpi
#export CC=mpicc CXX=mpicxx FC=mpifc
#export MPI_HOME=/lustre/sw/intel/compilers_and_libraries_2018.5.274/linux/mpi
#export MPI_LIB=$MPI_HOME/lib64
#export MPI_INCLUDE=$MPI_HOME/include64

# mpt
export CC=mpicc CXX=mpicxx FC=mpif90
export MPI_HOME=/opt/hpe/hpc/mpt/mpt-2.22
export MPI_INCLUDE=$MPI_HOME/include

# openmpi
#export CC=mpicc CXX=mpicxx FC=mpif90
#export MPI_HOME=/lustre/sw/openmpi/4.1.0-ucx-gcc8/
#export MPI_INCLUDE=$MPI_HOME/include

export LD_LIBRARY_PATH=/lustre/home/shared/dc025/binutils-devel/usr/lib64:$LD_LIBRARY_PATH

module use /lustre/home/shared/dc025/local_spack/modules/linux-rhel8-nehalem
module load libunwind-1.4.0-gcc-8.2.0-lwtjofw
module load libdwarf-20180129-gcc-8.2.0-scjswoi
module load libelf-0.8.13-gcc-8.2.0-ryuupys
module load papi-6.0.0.1-gcc-8.2.0-24pppx5
module load libxml2-2.9.10-gcc-8.2.0-4rlgln6
module load libiconv-1.16-gcc-8.2.0-kfblkpb
module load zlib-1.2.11-gcc-8.2.0-twb2vub

module load autotools
module load libtool

export UNWIND_HOME=/lustre/home/shared/dc025/local_spack/installs/linux-rhel8-nehalem/gcc-8.2.0/libunwind-1.4.0-lwtjofwhrfdft3dnizv4l7raaqga4ejy
export DWARF_HOME=/lustre/home/shared/dc025/local_spack/installs/linux-rhel8-nehalem/gcc-8.2.0/libdwarf-20180129-scjswoijrrikob5yxsnnwnnfvzueleh3
export ELF_HOME=/lustre/home/shared/dc025/local_spack/installs/linux-rhel8-nehalem/gcc-8.2.0/libelf-0.8.13-ryuupysfvtojdkfvkyzetqpp2fydqvpx
export PAPI_HOME=/lustre/home/shared/dc025/local_spack/installs/linux-rhel8-nehalem/gcc-8.2.0/papi-6.0.0.1-24pppx5vbsh23dchfsxsyfcrukdcgun5
export BINUTILS_HOME=/lustre/home/shared/dc025/binutils-devel/usr

mkdir -p $TARGET_EXTRAE
cd extrae-devel-3.8.4rc1
autoreconf -fiv
./configure --prefix=$TARGET_EXTRAE --with-mpi=$MPI_HOME --with-mpi-libs=$MPI_LIB --with-mpi-headers=$MPI_INCLUDE --with-unwind=$UNWIND_HOME --with-unwind-headers=$UNWIND_HOME/include --with-unwind-libs=$UNWIND_HOME/lib --without-dyninst --with-dwarf=$DWARF_HOME --with-elf=$ELF_HOME --with-papi=$PAPI_HOME --with-binutils=$BINUTILS_HOME
make
make install
