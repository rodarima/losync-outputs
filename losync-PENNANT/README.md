# Taskified PENNANT
This repository contains an unfinished version of the [PENNANT](https://github.com/lanl/PENNANT) software modified to use OmpSs-2 tasks and TAMPI.

## Approach
Each OpenMP loop in the ```Hydro::doCycle()``` function has been replaced with a loop generating a ```#pragma oss task```, one per "chunk".

The local sum and MPI operations in ```Mesh::sumToPoints()``` now use tasks. An MPI tag offset has been introduced to prevent MPI message matching errors between subsequent calls to ```Mesh::sumToPoints()```.

The call to ```Driver::calcGlobalDt()``` in the main loop is now a task. Note, this is a global synchronisation point for all tasks and processes as it contains a call to ```Parallel::globalMinLoc()```. This performs a collective ```MPI_Allreduce()``` to determine the minimum time step delta (```dt```) across all processes. All tasks within a process must have finished their respective ```dt``` calculations prior to this.

## Building
A sample build script is provided as ```do-build.sh```.

Preexisting installs of OmpSs-2 and TAMPI are required. The code has been tested with GCC 5.5.0. Newer GCC installations have been found to give build errors believed to be the result of C++11 incompatibilities with the Mercurium source-to-source compiler.

Environment variables ```TAMPI_HOME``` and ```GCC_HOME``` must point to the install directories of TAMPI and GCC respectively. Running ```make``` will then produce a build:

```
export TAMPI_HOME=/path/to/tampi
export GCC_HOME=/path/to/gcc-5.5.0
make
```

Optimisations and debug options can be enabled/disabled by modifying the ```Makefile``` directly as described in the PENNANT documentation.

## Running

Launch with ```mpirun``` or equivalent. Sample input files are available in ```test``` directory.

```
$ mpirun -n 2 ./build/pennant test/sedovsmall/sedovsmall.pnt
********************
Running PENNANT v0.9
********************

Running on 2 MPI PE(s)
--- Mesh Information ---
Points:  100
Zones:  81
Sides:  324
Edges:  189
Side chunks:  21
Point chunks:  8
Zone chunks:  6
Chunk size:  16
------------------------
Energy check:  total energy  =   2.467991e-01
(internal =   2.467991e-01, kinetic =   0.000000e+00)
End cycle      1, time = 2.50000e-03, dt = 2.50000e-03, wall = 1.11940e-02
dt limiter: Initial timestep
End cycle     10, time = 2.85593e-02, dt = 2.58849e-03, wall = 1.09181e-01
dt limiter: PE 0, Hydro dV/V limit for z = 0

Run complete
cycle =     10,         cstop =     10
time  =   2.855932e-02, tstop =   1.000000e+00

************************************
hydro cycle run time=   1.204391e-01
************************************
Energy check:  total energy  =   2.512181e-01
(internal =   1.874053e-01, kinetic =   6.381282e-02)
Writing .xy file...
```

## TODO
The current implementation contains several ```#pragma oss taskwait``` directives to synchronise operation. Without these ```taskwaits```, the code does not function (segfaults, incorrect outputs, etc.). It is suspected an error exists in the dependency declarations, as well as potentially further issues with MPI message matching.

It is unknown whether tasks will yield performance benefits given the ```Parallel::globalMinLoc()``` synchronisation point. A potential mitigation would be using a fixed ```dt``` across all time steps, rather than recalculating each cycle. Alternatively, a previous cycle's time step size could be used to estimate sizes for subsequent time steps.
