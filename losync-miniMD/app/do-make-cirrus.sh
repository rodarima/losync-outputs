#!/bin/bash
set -e

module load gcc/8.2.0
#export PREFIX=/lustre/home/shared/dc025/extrae-3.7.1-test
#export PREFIX=/lustre/home/shared/dc025/ompss2-2020.11.1-openmpi
export PREFIX=/lustre/home/shared/dc025/ompss2-2020.11.1-mpt
export TAMPI_HOME=${PREFIX}/tampi-install
export NANOS6_HOME=${PREFIX}/nanos6-install
export MCXX_HOME=${PREFIX}/mercurium
export PATH=${MCXX_HOME}/bin:$PATH

# Ompss-2
#CFLAGS="--ompss-2 -I${TAMPI_HOME}/include -O3 -DMPICH_IGNORE_CXX_SEEK -DNOCHUNK -mavx  -UUSE_SIMD -DUSE_TAMPI -DUSE_TASKS -DPRECISION=2"
# Ompss-2 Traced
CFLAGS="--ompss-2 -I${TAMPI_HOME}/include -O3 -g -DMPICH_IGNORE_CXX_SEEK -DNOCHUNK -mavx  -UUSE_SIMD -DUSE_TAMPI -DUSE_TASKS -DPRECISION=2"
# Ompss-2 debug
#CFLAGS="--ompss-2 -I${TAMPI_HOME}/include -O0 -g -fsanitize=address -DMPICH_IGNORE_CXX_SEEK -DNOCHUNK -mavx -UUSE_SIMD -DUSE_TAMPI -DUSE_TASKS -DPRECISION=2"
# No Ompss-2/TAMPI
#CFLAGS="-O3 -DMPICH_IGNORE_CXX_SEEK -DNOCHUNK -mavx  -UUSE_SIMD -UUSE_TAMPI -UUSE_TASKS -DPRECISION=2"

SOURCE_FILES="ljs.cpp input.cpp integrate.cpp atom.cpp force_lj.cpp force_eam.cpp neighbor.cpp thermo.cpp comm.cpp timer.cpp output.cpp setup.cpp"

# intel mpi
#module load intel-mpi-18
#cmd="mpigxx -cxx=mcxx ${CFLAGS} ${SOURCE_FILES} -o miniMD_boxes -ltampi -L${TAMPI_HOME}/lib"

# openmpi
#module load openmpi/4.1.0
#export OMPI_MPICC=mcc
#export OMPI_MPICXX=mcxx
#cmd="mpicxx ${CFLAGS} ${SOURCE_FILES} -o miniMD_boxes -lmpi_mpifh -ltampi -L${TAMPI_HOME}/lib"

# mpt
module load mpt
export MPICC_CC=mcc
export MPICXX_CXX=mcxx
cmd="mpicxx ${CFLAGS} ${SOURCE_FILES} -o miniMD_boxes -ltampi -L${TAMPI_HOME}/lib"

rm -f *.o
rm -f miniMD_boxes
echo $cmd
eval $cmd
