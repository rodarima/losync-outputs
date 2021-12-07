#!/bin/bash
set -e
export MERCURIUM=/lustre/home/shared/dc025/ompss2-2021.06-mpt+extrae/mercurium
export PKG_CONFIG_PATH=/home/shared/dc025/prog/sqlite3/lib/pkgconfig
module load gcc/8.2.0
module load libtool

module use /lustre/home/shared/dc025/local_spack/modules/linux-rhel8-nehalem
module load spack
module load gperf-3.1-gcc-8.2.0-7q3dayj
module load bison-3.4.2-gcc-8.2.0-saivm44
module load flex-2.6.4-gcc-8.2.0-zlwjqca

mkdir -p $MERCURIUM
cd ompss-2-releases/mcxx
#autoreconf -fiv
#make clean
#./configure --prefix=$MERCURIUM --enable-ompss-2 --with-nanos6=/lustre/home/shared/dc025/ompss2-2021.06-mpt+extrae/nanos6-install
./configure --prefix=$MERCURIUM --enable-ompss-2 --enable-nanos6-bootstrap
make
make install
