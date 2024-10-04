#!/bin/bash

#eval `../../../../../cime/CIME/Tools/get_case_env`
E3SM_ROOT="/pscratch/sd/z/zhan391/DARPA_project/e3sm_model/code/E3SM"

${E3SM_ROOT}/cime/tools/configure && source .env_mach_specific.sh

export FC=gfortran
export LIB_NETCDF=${NETCDF_DIR}/lib
export INC_NETCDF=${NETCDF_DIR}/include
make
