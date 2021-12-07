#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=miniMD
#SBATCH --time=0:10:0
#SBATCH --exclusive
#SBATCH --nodes=28
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=18

# Replace [budget code] below with your budget code (e.g. t01)
#SBATCH --account=dc025
# We use the "standard" partition as we are running on CPU nodes
#SBATCH --partition=standard
# We use the "standard" QoS as our runtime is less than 4 days
#SBATCH --qos=standard

# Nanos6 >=2022.11.1 loader fails if papi module is not loaded
module use /lustre/home/shared/dc025/local_spack/modules/linux-rhel8-nehalem
module load spack
module load papi-6.0.0.1-gcc-8.2.0-24pppx5

#module load intel-mpi-18
#module load openmpi/4.1.0
module load mpt

module load gcc/8.2.0
module load gdb
#module load forge
#export LD_LIBRARY_PATH=/lustre/sw/allinea/Forge-20/lib/64:${LD_LIBRARY_PATH}
#export PREFIX=/lustre/home/shared/dc025/extrae-3.7.1-test
#export PREFIX=/lustre/home/shared/dc025/ompss2-2020.11.1
#export PREFIX=/lustre/home/shared/dc025/ompss2-2020.11.1-openmpi
#export PREFIX=/lustre/home/shared/dc025/ompss2-2021.06-mpt
export PREFIX=/lustre/home/shared/dc025/ompss2-2021.06-mpt+extrae
export TAMPI_HOME=${PREFIX}/tampi-install
export NANOS6_HOME=${PREFIX}/nanos6-install
export EXTRAE_HOME=${PREFIX}/extrae
export LD_LIBRARY_PATH=${TAMPI_HOME}/lib:${NANOS6_HOME}/lib:${EXTRAE_HOME}/lib:$LD_LIBRARY_PATH

# Needed for generating ctf traces
module load python-3.7.7-gcc-8.2.0-ndlvqr4
export PATH=${PREFIX}/nanos6-install/bin:/lustre/home/shared/dc025/babeltrace2-2.0.3/bin:$PATH
export PYTHONPATH=/lustre/home/shared/dc025/babeltrace2-2.0.3/lib/python3.7/site-packages:${PREFIX}/nanos6-install/share/doc/nanos6/scripts/ctf/plugins:$PYTHONPATH

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=1

#export NANOS6=extrae
#export EXTRAE_ON=1
export NANOS6_CONFIG=$PWD/nanos6.toml

# Confirm which MPI we're using
#ldd ./miniMD_boxes

#srun --cpu-bind=verbose --cpu-bind=map_cpu:$(seq -s, 0 36) ./miniMD_boxes --half_neigh 0 -i in.losync.cirrus
#srun --cpu-bind=mask_cpu:0x3FFFF,0xFFFFC0000 ./miniMD_boxes --half_neigh 0 -i in.losync.cirrus
srun --cpu-bind=threads ./miniMD_boxes --half_neigh 0 -i in.losync.cirrus
#srun --cpu-bind=threads ./trace.sh ./miniMD_boxes --half_neigh 0 -i in.losync.cirrus
#srun ./miniMD_boxes --half_neigh 0 -i in.losync.cirrus

