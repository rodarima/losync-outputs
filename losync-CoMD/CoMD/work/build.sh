#!/bin/bash

module load mpt
module swap gcc/6.3.0 gcc/8.2.0
export NANOS6_HOME=/lustre/home/shared/dc025/ompss2-2021.06-mpt
export MCXX_HOME=/lustre/home/shared/dc025/ompss2-2021.06-mpt/mercurium
export PATH=${MCXX_HOME}/bin:$PATH
export MPICC_CC=mcc

cd ../src-tampi
make clean
if [ $1 = "tampi" ]; then
    make tampi
else
    make
fi
