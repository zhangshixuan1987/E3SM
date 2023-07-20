# This file is for user convenience only and is not used by the model
# Changes to this file will be ignored and overwritten
# Changes to the environment should be made in env_mach_specific.xml
# Run ./case.setup --reset to regenerate this file
. /usr/share/lmod/8.3.1/init/sh
module unload cray-hdf5-parallel cray-netcdf-hdf5parallel cray-parallel-netcdf PrgEnv-gnu PrgEnv-intel PrgEnv-nvidia PrgEnv-cray PrgEnv-aocc intel intel-oneapi cudatoolkit craype-accel-nvidia80 craype-accel-host perftools-base perftools darshan
module load PrgEnv-gnu/8.3.3 gcc/11.2.0 cray-libsci/23.02.1.1 craype-accel-host craype/2.7.19 cray-mpich/8.1.24 cray-hdf5-parallel/1.12.2.3 cray-netcdf-hdf5parallel/4.9.0.3 cray-parallel-netcdf/1.12.3.3 cmake/3.24.3
export MPICH_ENV_DISPLAY=1
export MPICH_VERSION_DISPLAY=1
export OMP_STACKSIZE=128M
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export HDF5_USE_FILE_LOCKING=FALSE
export PERL5LIB=/global/cfs/cdirs/e3sm/perl/lib/perl5-only-switch
export FI_CXI_RX_MATCH_MODE=software
export MPICH_COLL_SYNC=MPI_Bcast
export ADIOS2_DIR=/global/cfs/cdirs/e3sm/3rdparty/adios2/2.8.3.patch/cray-mpich-8.1.15/gcc-11.2.0
export LD_LIBRARY_PATH=/pscratch/sd/z/zhan391/DARPA_project/e3sm_model/code/deepONet/cpython/lib/:/pscratch/sd/z/zhan391/DARPA_project/e3sm_model/code/deepONet/cpython/lib/:/opt/cray/pe/gcc/11.2.0/snos/lib64:/opt/cray/libfabric/1.15.2.0/lib64::
export PATH=/pscratch/sd/z/zhan391/DARPA_project/e3sm_model/code/deepONet/cpython/bin:/pscratch/sd/z/zhan391/DARPA_project/e3sm_model/code/deepONet/cpython/bin:/global/common/software/nersc/pm-2022q4/spack/linux-sles15-zen/cmake-3.24.3-k5msymx/bin:/opt/cray/pe/parallel-netcdf/1.12.3.3/bin:/opt/cray/pe/netcdf-hdf5parallel/4.9.0.3/bin:/opt/cray/pe/hdf5-parallel/1.12.2.3/bin:/opt/cray/pe/hdf5/1.12.2.3/bin:/opt/cray/pe/mpich/8.1.24/ofi/gnu/9.1/bin:/opt/cray/pe/mpich/8.1.24/bin:/opt/cray/pe/craype/2.7.19/bin:/opt/cray/pe/gcc/11.2.0/bin:/usr/common/software/python/3.8-anaconda-2020.11/bin:/global/common/software/e3sm/anaconda_envs/base/bin:/global/common/software/e3sm/anaconda_envs/base/condabin:/global/common/software/nersc/bin:/opt/cray/libfabric/1.15.2.0/bin:/global/homes/z/zhan391/.local/bin:/usr/local/bin:/usr/bin:/bin:/usr/lib/mit/bin:/opt/cray/pe/bin