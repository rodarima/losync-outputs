#!/bin/bash
module use /lustre/home/shared/dc025/local_spack/modules/linux-rhel8-nehalem
module load spack
module load python-3.7.7-gcc-8.2.0-ndlvqr4
export PATH=/lustre/home/shared/dc025/ompss2-2021.06-mpt/nanos6-install/bin:/lustre/home/shared/dc025/babeltrace2-2.0.3/bin:$PATH
export PYTHONPATH=/lustre/home/shared/dc025/babeltrace2-2.0.3/lib/python3.7/site-packages:/lustre/home/shared/dc025/ompss2-2021.06-mpt/nanos6-install/share/doc/nanos6/scripts/ctf/plugins:$PYTHONPATH
