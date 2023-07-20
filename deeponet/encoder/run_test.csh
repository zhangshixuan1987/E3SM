#!/bin/csh

#source .env_mach_specific.csh
eval `../../cime/CIME/Tools/get_case_env`

set topdir = "/pscratch/sd/z/zhan391/DARPA_project/e3sm_model/code/E3SM"
setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:${topdir}/pytorch-lib/lib/python3.10/site-packages/torch/lib:${topdir}/pytorch-fortran/gnu/install/lib"
setenv PYTHONPATH "${topdir}/pytorch-lib/lib/python3.10/site-packages"


