# Taskified CoMD
This repository contains a version of the [CoMD](https://github.com/ECP-copa/CoMD) software modified to use OmpSs-2 tasks and TAMPI.

## Approach
The taskified code is located in ```CoMD/src-tampi```. The communication pattern
has been changed so that each processor communicates directly with each of its
neighbours. The original code only communicated with its orthogonally adjacent
neighbours, and dealt with corner halo cells by passing them through multiple
communications. The calculation and communication phases are now processed in a
column-by-column fashion instead of box-by-box. This is because the work
involved in a single box was too small for tasks to be useful

## Building on the Cirrus HPC Service
A sample build script is provided in ```CoMD/work/build.sh```.

Preexisting installs of OmpSs-2 and TAMPI are required. 

## Running on the Cirrus HPC Service
A sample job script is provided in ```CoMD/work/submit.slurm```.

A combined script to build and run the code is provided in ```CoMD/work/run.sh```.
This script also processes the output files and convert relevant information
into a more readable format.

Run-time parameters can be found in ```CoMD/work/parameters.config```.
