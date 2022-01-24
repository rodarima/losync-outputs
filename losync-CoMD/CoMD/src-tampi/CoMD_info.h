#ifndef CoMD_info_hpp
#define CoMD_info_hpp

#define CoMD_VARIANT "CoMD-openmp-mpi"
#define CoMD_HOSTNAME "cirrus-login2"
#define CoMD_KERNEL_NAME "'Linux'"
#define CoMD_KERNEL_RELEASE "'4.18.0-147.8.1.el8_1.x86_64'"
#define CoMD_PROCESSOR "'x86_64'"

#define CoMD_COMPILER "'/opt/hpe/hpc/mpt/mpt-2.22/bin/mpicc'"
#define CoMD_COMPILER_VERSION "'mcxx 2.3.0 (unknown revision)'"
#define CoMD_CFLAGS "'-std=c99 --ompss-2 -DDOUBLE -DDO_MPI   -I/lustre/home/shared/dc025/ompss2-2021.06-mpt/tampi-install/include/ '"
#define CoMD_LDFLAGS "'-ltampi-c -L/lustre/home/shared/dc025/ompss2-2021.06-mpt/tampi-install/lib/ -lm '"

#endif
