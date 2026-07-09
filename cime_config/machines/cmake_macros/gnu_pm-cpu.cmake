string(APPEND CONFIG_ARGS " --host=cray")

if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS
    " -DHAVE_NANOTIME"
    " -DBIT64"
    " -DHAVE_SLASHPROC"
    " -DHAVE_GETTIMEOFDAY"
  )
endif()

string(APPEND SLIBS
  " -L$ENV{CRAY_HDF5_PARALLEL_PREFIX}/lib"
  " -lhdf5_hl"
  " -lhdf5"
  " -L$ENV{CRAY_NETCDF_HDF5PARALLEL_PREFIX}/lib"
  " -L$ENV{CRAY_PARALLEL_NETCDF_PREFIX}/lib"
  " -lpnetcdf"
  " -lnetcdf"
  " -lnetcdff"
)

# Use Intel oneMKL instead of the system BLAS/LAPACK libraries.
# The system /usr/lib64/libblas.so depends on libgfortran.so.4,
# while the GNU 12 and PyTorch libraries use libgfortran.so.5.
string(APPEND SLIBS
  " -L/opt/intel/oneapi/mkl/latest/lib/intel64"
  " -lmkl_rt"
  " -lpthread"
  " -lm"
  " -ldl"
)

set(CXX_LINKER "FORTRAN")

set(NETCDF_PATH "$ENV{CRAY_NETCDF_HDF5PARALLEL_PREFIX}")
set(NETCDF_C_PATH "$ENV{CRAY_NETCDF_HDF5PARALLEL_PREFIX}")
set(NETCDF_FORTRAN_PATH "$ENV{CRAY_NETCDF_HDF5PARALLEL_PREFIX}")
set(HDF5_PATH "$ENV{CRAY_HDF5_PARALLEL_PREFIX}")
set(PNETCDF_PATH "$ENV{CRAY_PARALLEL_NETCDF_PREFIX}")

if (NOT DEBUG)
  string(APPEND CFLAGS
    " -O2"
    " -g"
    " -fPIC"
  )

  string(APPEND FFLAGS
    " -O2"
    " -g"
    " -fbounds-check"
    " -fPIC"
  )
endif()

string(APPEND CXX_LIBS " -lstdc++")
string(APPEND CXXFLAGS " -D_GLIBCXX_USE_CXX11_ABI=1")

# PyTorch-Fortran module/include paths.
string(APPEND FFLAGS
  " -I$ENV{TORCH_DIR}/torch/include"
  " -I$ENV{TORCH_DIR}/torch/include/mod_files"
)

# PyTorch-Fortran libraries and Intel ITT/JIT symbols required by
# libtorch_cpu.so.
string(APPEND SLIBS
  " -L$ENV{TORCH_DIR}/torch/lib"
  " -Wl,--no-as-needed"
  " -lpytorch_fort_proxy"
  " -lpytorch_proxy"
  " -Wl,--whole-archive"
  " $ENV{TORCH_DIR}/build/lib/libittnotify.a"
  " -Wl,--no-whole-archive"
  " -Wl,--as-needed"
)

set(MPICC "cc")
set(MPICXX "CC")
set(MPIFC "ftn")

set(SCC "gcc")
set(SCXX "g++")
set(SFC "gfortran")
