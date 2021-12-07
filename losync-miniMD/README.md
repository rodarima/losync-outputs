# Taskified miniMD
This repository contains a version of the [Mantevo miniMD](https://github.com/Mantevo/miniMD) software modified to use OmpSs-2 tasks and TAMPI.

## Approach
A high-level overview is given in ```doc/minimd-approach-overview.pdf```. Further to this, the ```communicate```, ```exchange``` and ```borders``` MPI steps have been decomposed into individual ```pack```, ```send```, ```recv``` and ```unpack``` tasks, with the ```pack+send``` and ```recv+unpack``` tasks being merged in cases where maintaining them seperately afforded no additional parallelism.

Several implementations of ```Comm::communicate()```, ```Comm::exchange()``` and ```Comm::borders()``` exist and can be switched between by uncommmenting the relevant lines in each function. For example:

```
//communicate_nonblocking_alltasks_tampi_iwait(atoms[box_index], box_index);
communicate_nonblocking_neighbourtasks_tampi_iwaitall(atoms[box_index], box_index);
```

to:

```
communicate_nonblocking_alltasks_tampi_iwait(atoms[box_index], box_index);
//communicate_nonblocking_neighbourtasks_tampi_iwaitall(atoms[box_index], box_index);
```

switches an implementation of the ```communicate()``` function that uses a single task per process neighbour (i.e. 3 non-blocking MPI messages per task) and the ```TAMPI_Iwaitall()``` function with one that uses a single task per MPI message (i.e. 3 tasks per neighbour) and the ```TAMPI_Iwait()``` function. Additional implementations are documented within the code.

## Building on the Cirrus HPC Service
A sample build script is provided as ```do-make-cirrus.sh```.

Preexisting installs of OmpSs-2 and TAMPI are required. The ```mpt``` MPI implementation is recommended as it performs significantly better in miniMD benchmarks, however builds can be performed against ```intel-mpi``` if required. Building against ```openmpi``` is untested.

To facilitate building against OmpSs-2 and TAMPI, the provided Makefiles are not used. Compiler wrapper ```mpicxx``` is used directly:

```
CFLAGS="--ompss-2 -I${TAMPI_HOME}/include -O3 -g -DMPICH_IGNORE_CXX_SEEK -DNOCHUNK -mavx  -UUSE_SIMD -DUSE_TAMPI -DUSE_TASKS -DPRECISION=2"

SOURCE_FILES="ljs.cpp input.cpp integrate.cpp atom.cpp force_lj.cpp force_eam.cpp neighbor.cpp thermo.cpp comm.cpp timer.cpp output.cpp setup.cpp"

module load mpt
export MPICC_CC=mcc
export MPICXX_CXX=mcxx
mpicxx ${CFLAGS} ${SOURCE_FILES} -o miniMD_boxes -ltampi -L${TAMPI_HOME}/lib
```

## Running on the Cirrus HPC Service
A sample job script is provided as ```cirrus_jobscript.slurm```.

Input file ```in.losync.cirrus``` includes the standard benchmarking case used during development. Settings are documented within the file.

The "boxes per process" setting adjusts the number of tasks produced by each process. Default is 17: each Cirrus CPU contains 18 cores which allows for 1 box per core + 1 core for the task generation thread.

The "processes in x-direction" and "processes in z-direction" settings control the process grid and must be adjusted to match the process count for each job. For example, for a ```--nodes=28``` and ```--tasks-per-node=2``` job, the total number of processes will be 28x2 = 56. A suitable "processes in x-direction" and "processes in z-direction" would therefore be 7 and 8 respectively (7x8 = 56). There is no "processes in y-direction" as "boxes per process" is equivalent.

When launching, input parameter ```--half_neigh 0``` should be set as the half neighbour mode from the original miniMD is not supported by the taskified version. Use of Hyperthreading is not recommended and can be disabled with the ```--cpu-bind=threads``` option to ```srun```:

```srun --cpu-bind=threads ./miniMD_boxes --half_neigh 0 -i in.losync.cirrus```

Task runtime settings are controlled by the ```nanos6.toml``` file. Given the high number of tasks generated, setting ```policy = "busy"``` under ```[cpumanager]``` has been found to give marginal performance improvements.

## Performance
Run times of the taskified and original, unmodified code on Cirrus are provided in ```doc/minimd_cirrus_runtimes.ods```.

Performance can be measured by examining the ```# Performance Summary:``` lines of the output:

```
# Starting dynamics ...
# Timestep T U P Time
0 1.440000e+00 -6.773368e+00 -5.019671e+00  0.000
100 8.071170e-01 -5.831598e+00 -1.579241e-01  4.543


# Performance Summary:
# MPI_proc OMP_threads nsteps natoms t_total t_force t_neigh t_comm t_other performance perf/thread grep_string t_extra
56 1 100 907924 4.542873 1.463951 0.286836 2.782119 0.009967 19985680.080783 356887.144300 PERF_SUMMARY 0.388933
```

The ```t_total``` value gives the run time in seconds of the main loop, here 4.542873s.

## TODO
The best performing implementation of the taskified version (```communicate_nonblocking_neighbourtasks_tampi_iwaitall()```) completes a 61x61x61 problem with 2500 time steps in approximately twice the runtime of unmodified miniMD. Profiling suggests the source of the additional runtime are the ```communicate_send``` and ```communicate_recv``` tasks but investigations are on-going.

A single, best approach for ```Comm::communicate()```, ```Comm::exchange()``` and ```Comm::borders()``` should be decided and the redundant implementations removed from the code. If different approaches are optimal for different cases, it should be possible to select an implementation from the input file.
