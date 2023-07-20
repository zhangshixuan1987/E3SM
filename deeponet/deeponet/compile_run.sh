#/bin/bash

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/users/kshukla1/WORK_RAJ/E3SM_INTEG/pytorch-fortran/gnu/install/lib
gfortran resnet_forward.f90 -o exec -g -fbounds-check -I/users/kshukla1/WORK_RAJ/E3SM_INTEG/pytorch-fortran/gnu/install/include/mod_files -L/users/kshukla1/WORK_RAJ/E3SM_INTEG/pytorch-fortran/gnu/install/lib -lpytorch_fort_proxy -lpytorch_proxy

