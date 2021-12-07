# Building TAMPI and OmpSs-2 on Cirrus HPC Service
Steps for building on the [Cirrus HPC Service](https://www.cirrus.ac.uk).

## Prerequisites
Using the [Spack](https://spack.readthedocs.io) package manager and appropriate version of GCC:

```
module load spack
module load gcc/8.2.0
```

install the following:

```
autoconf-2.69-gcc-8.2.0-bkc32sr
automake-1.16.2-gcc-8.2.0-lyeujjm
binutils-2.34-gcc-8.2.0-ezr3mdx
bison-3.4.2-gcc-8.2.0-saivm44
bzip2-1.0.8-gcc-8.2.0-zagvcfz
diffutils-3.7-gcc-8.2.0-hcrayh4
elfutils-0.179-gcc-8.2.0-msoboez
expat-2.2.9-gcc-8.2.0-3le7aqk
findutils-4.6.0-gcc-8.2.0-rvzzefr
flex-2.6.4-gcc-8.2.0-zlwjqca
gdb-9.1-gcc-8.2.0-kedc2pr
gdbm-1.18.1-gcc-8.2.0-b6x5e4d
gettext-0.20.2-gcc-8.2.0-upr3dit
gperf-3.1-gcc-8.2.0-7q3dayj
graphviz-2.42.2-gcc-8.2.0-brphsrb
help2man-1.47.11-gcc-8.2.0-hhgtxcx
hwloc-2.2.0-gcc-8.2.0-vosjebz
libbsd-0.10.0-gcc-8.2.0-nsxm6bt
libdwarf-20180129-gcc-8.2.0-scjswoi
libelf-0.8.13-gcc-8.2.0-ryuupys
libffi-3.3-gcc-8.2.0-32vqjac
libiconv-1.16-gcc-8.2.0-kfblkpb
libpciaccess-0.13.5-gcc-8.2.0-vdhcxnn
libsigsegv-2.12-gcc-8.2.0-wpybls2
libtool-2.4.6-gcc-8.2.0-3l4qrtz
libunwind-1.4.0-gcc-8.2.0-lwtjofw
libxml2-2.9.10-gcc-8.2.0-4rlgln6
m4-1.4.18-gcc-8.2.0-lxrtr33
ncurses-6.2-gcc-8.2.0-y5qwbjv
numactl-2.0.12-gcc-8.2.0-pl5uzim
openssl-1.1.1g-gcc-8.2.0-oj2lnyu
papi-6.0.0.1-gcc-8.2.0-24pppx5
parallel-20190222-gcc-8.2.0-nqmjdbj
pcre-8.44-gcc-8.2.0-abbeas3
perl-5.30.2-gcc-8.2.0-sl4nxwj
pkgconf-1.7.3-gcc-8.2.0-okmsxq7
python-3.7.7-gcc-8.2.0-ndlvqr4
readline-8.0-gcc-8.2.0-eisajws
sqlite-3.31.1-gcc-8.2.0-l4pzb7o
swig-4.0.1-gcc-8.2.0-xncdjoc
tar-1.32-gcc-8.2.0-h3ysjua
texinfo-6.5-gcc-8.2.0-o3lz3i5
util-macros-1.19.1-gcc-8.2.0-s7fgk5m
xz-5.2.5-gcc-8.2.0-gn3o3oj
zlib-1.2.11-gcc-8.2.0-twb2vub
```

In order to perform local installs with Spack, the ```~/.spack/config.yaml``` file must be updated/created with ```install_tree``` and all other relevant paths set to a directory writeable by the current user. A sample ```config.yaml``` file is provided in this repository.

## Acquiring Source Code

TAMPI: https://github.com/bsc-pm/tampi

OmpSs-2: https://github.com/bsc-pm/ompss-2-releases

Extrae (optional, used for tracing): https://github.com/bsc-performance-tools/extrae

Babeltrace2 (optional, used for tracing): https://babeltrace.org/#bt2-get

## Guide to Build Scripts

Sample ```build-*.sh``` scripts are provided. Relevant paths will need updated before use. The ```mpt``` MPI implementation is recommended as it performs significantly better in benchmarks of the taskified miniMD mini-app, however builds can be performed against ```intel-mpi``` and ```openmpi``` by uncommenting the relevant lines.

If required, Babeltrace2 and Extrae must be built first.

TAMPI can be built in any order.

Mercurium depends on Nanos6 and Nanos6 depends on Mercurium. To break this "chicken and egg" dependence:

* Build Mercurium with the ```--enable-nanos6-bootstrap``` flag
* Build Nanos6 against this Mercurium install
* Rebuild Mercurium against this Nanos6 install
* Finally, rebuild Nanos6 against the Mercurium rebuild

Nanos6 builds give errors with ```make -j <num>``` so must be performed serially. On Cirrus, this should be expected to take several hours.

