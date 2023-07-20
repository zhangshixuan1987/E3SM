#!/bin/bash -l 
source .env_mach_specific.sh  
topdir=/pscratch/sd/z/zhan391/DARPA_project/e3sm_model/code/E3SM/pytorch-lib
export PATH="${PATH}:${topdir}/lib/python3.10/site-packages"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${topdir}/lib/python3.10/site-packages/torch/lib"
export PYTHONPATH="${PYTHONPATH}:${topdir}/lib/python3.10/site-packages"
echo $PYTHONPATH

