# This file is for user convenience only and is not used by the model
# Changes to this file will be ignored and overwritten
# Changes to the environment should be made in env_mach_specific.xml
# Run ./case.setup --reset to regenerate this file

source /global/common/software/e3sm/anaconda_envs/base/etc/profile.d/conda.sh
source /global/common/software/e3sm/anaconda_envs/base/etc/profile.d/mamba.sh

export E3SMU_SCRIPT="/global/common/software/e3sm/anaconda_envs/load_e3sm_unified_1.8.1_pm-cpu.sh"
export E3SMU_MACHINE="pm-cpu"

export E3SMU_MPI="SYSTEM"
mamba activate e3sm_unified_1.8.1_pm-cpu

source /global/common/software/e3sm/anaconda_envs/spack/spack_for_mache_1.10.0/share/spack/setup-env.sh
spack env activate e3sm_unified_1_8_1_pm-cpu_gnu_mpich
module rm PrgEnv-gnu
module rm PrgEnv-nvidia
module rm cudatoolkit
module rm craype-accel-nvidia80
module rm craype-accel-host
module rm perftools-base
module rm perftools
module rm darshan

module load PrgEnv-gnu/8.3.3
module load gcc/11.2.0
module load craype-accel-host

module load craype
module rm cray-mpich
module load libfabric/1.15.2.0
module load cray-mpich/8.1.22

module load cmake/3.22.0

export MPICH_ENV_DISPLAY=1
export MPICH_VERSION_DISPLAY=1
export OMP_STACKSIZE=128M
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export HDF5_USE_FILE_LOCKING=FALSE
export PERL5LIB=/global/cfs/cdirs/e3sm/perl/lib/perl5-only-switch
export FI_CXI_RX_MATCH_MODE=software

export HDF5_USE_FILE_LOCKING=FALSE
