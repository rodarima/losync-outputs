#!/bin/bash
set -e
module use /lustre/home/shared/dc025/local_spack/modules/linux-rhel8-nehalem
module load spack
module load elfutils-0.179-gcc-8.2.0-msoboez
module load swig-4.0.1-gcc-8.2.0-xncdjoc
module load python-3.7.7-gcc-8.2.0-ndlvqr4
module load gettext-0.20.2-gcc-8.2.0-upr3dit
module load flex-2.6.4-gcc-8.2.0-zlwjqca
module load bison-3.4.2-gcc-8.2.0-saivm44

module load gcc/8.2.0

cd babeltrace2-2.0.3
./configure --prefix=/lustre/home/shared/dc025/babeltrace2-2.0.3 --enable-python-bindings --enable-python-plugins
make -j8
make install
