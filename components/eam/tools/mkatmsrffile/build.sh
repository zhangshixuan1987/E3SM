#!/bin/bash
eval `/pscratch/sd/z/zhan391/DARPA_project/e3sm_model/code/E3SM/cime/CIME/Tools/get_case_env`
export FC=gfortran
export NETCDF_ROOT=$NETCDF_DIR
make

