#!/bin/bash

set -ex

module load anaconda/python3

prefix=""

if [ $# -eq 0 ]; then
    echo 'Please enter walltime'
    exit 1
fi

if [ $2 = "tampi" ]; then
    ./build.sh tampi
    prefix="tampi_"
else
    ./build.sh
fi

timestamp=$(date +%Y-%m-%d_%H-%M-%S)

job="submit.slurm"
jobID_full=$(sbatch --time=$1 $job)
jobID=$(echo "$jobID_full" | tr -dc '0-9')

numDone=0

while [ $numDone -lt 1 ]; do
    if squeue --job $jobID | grep -q $jobID; then
        clear
        squeue --job $jobID
        sleep 10
    else
        mkdir yamls/"$timestamp"
        mv *.yaml yamls/"$timestamp"

        python output_parser.py slurm-"$jobID".out timing_"$timestamp".dat
        mv timing_"$timestamp".dat results/"$prefix""$timestamp".dat
        mv slurm-"$jobID".out results/"$prefix"CoMD_"$timestamp".dat
        echo "results/""$prefix""CoMD_""$timestamp"".dat"
        numDone=1
    fi
done
