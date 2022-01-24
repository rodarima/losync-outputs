#!/bin/bash

#SBATCH --job-name=CoMD
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=18

#SBATCH --account=dc025
#SBATCH --partition=standard
#SBATCH --qos=standard

module load mpt
module use /lustre/home/shared/dc025/local_spack/modules/linux-rhel8-nehalem
module load spack
module load papi-6.0.0.1-gcc-8.2.0-24pppx5

module swap gcc/6.3.0 gcc/8.2.0

export PREFIX=/lustre/home/shared/dc025/ompss2-2021.06-mpt
export TAMPI_HOME=${PREFIX}/tampi-install
export NANOS6_HOME=${PREFIX}/nanos6-install
export LD_LIBRARY_PATH=${TAMPI_HOME}/lib:${NANOS6_HOME}/lib:$LD_LIBRARY_PATH

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=18

cd $SLURM_SUBMIT_DIR

source parameters.config

# Launch the parallel job
srun  ../bin/CoMD-openmp-mpi -i$xproc -j$yproc -k$zproc -x $x -y $y -z $z -N $steps -n $printrate
