#!/bin/bash
set -e

module load gcc/8.2.0
module load mpt
#module load intel-mpi-18
#module load openmpi/4.1.0

CFLAGS="-O0 -g -DMPICH_IGNORE_CXX_SEEK -DNOCHUNK -mavx  -UUSE_SIMD -UUSE_TAMPI -UUSE_TASKS -DPRECISION=2 -Wl,--allow-multiple-definition"
#LFLAGS="-ldmallocthcxx"
#CFLAGS="-O0 -g -fsanitize=address -DMPICH_IGNORE_CXX_SEEK -DNOCHUNK -mavx  -UUSE_SIMD -UUSE_TAMPI -UUSE_TASKS -DPRECISION=2"
#CFLAGS="-O3 -DMPICH_IGNORE_CXX_SEEK -DNOCHUNK -mavx  -UUSE_SIMD -UUSE_TAMPI -UUSE_TASKS -DPRECISION=2"
SOURCE_FILES="ljs.cpp input.cpp integrate.cpp atom.cpp force_lj.cpp force_eam.cpp neighbor.cpp thermo.cpp comm.cpp timer.cpp output.cpp setup.cpp"

rm -f *.o
rm -f miniMD_boxes
cmd="mpicxx ${CFLAGS} ${SOURCE_FILES} -o miniMD_boxes ${LFLAGS}"
echo $cmd
eval $cmd
