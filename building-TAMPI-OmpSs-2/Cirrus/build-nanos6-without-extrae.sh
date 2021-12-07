#!/bin/bash
set -e
module use /lustre/home/shared/dc025/local_spack/modules/linux-rhel8-nehalem
module load spack
module load hwloc-2.2.0-gcc-8.2.0-vosjebz
module load libxml2-2.9.10-gcc-8.2.0-4rlgln6
module load libpciaccess-0.13.5-gcc-8.2.0-vdhcxnn
module load numactl-2.0.12-gcc-8.2.0-pl5uzim
module load boost/1.73.0
module load libunwind-1.4.0-gcc-8.2.0-lwtjofw
module load papi-6.0.0.1-gcc-8.2.0-24pppx5
module load elfutils-0.179-gcc-8.2.0-msoboez
module load graphviz-2.42.2-gcc-8.2.0-brphsrb
module load parallel-20190222-gcc-8.2.0-nqmjdbj
module load zlib-1.2.11-gcc-8.2.0-twb2vub
module load pkgconf-1.7.3-gcc-8.2.0-okmsxq7

module load autotools
module load libtool

module load mpt
#module load openmpi/4.1.0
#module load intel-20.4/mpi
module load gcc/8.2.0

TARGET=/lustre/home/shared/dc025/ompss2-2021.06-mpt/nanos6-install
MCXX_HOME=/lustre/home/shared/dc025/ompss2-2021.06-mpt/mercurium
mkdir -p $TARGET
cd ompss-2-releases/nanos6
autoreconf -fiv
#make clean
./configure --prefix=$TARGET --with-nanos6-mercurium=$MCXX_HOME --with-babeltrace2=/lustre/home/shared/dc025/babeltrace2-2.0.3
make all check
#make check
make install
