################################################################################
# Initial cache list for Cray systems with GPU (Piz Daint, Titan)
# Usage: cmake -C path/to/Cray_GPU.cmake ...
################################################################################

# The C++ compile wrapper CC already includes and links to these libraries.
# No need to look for them.
set(DCA_HDF5_IMPLICIT   TRUE CACHE INTERNAL "")
set(DCA_LAPACK_IMPLICIT TRUE CACHE INTERNAL "")
set(DCA_FFTW_IMPLICIT   TRUE CACHE INTERNAL "")

# MPIEXEC stuff for executing parallel tests.
set(MPIEXEC "aprun"
  CACHE FILEPATH "Executable for running MPI programs.")
set(MPIEXEC_NUMPROC_FLAG "-n"
  CACHE FILEPATH
  "Flag used by MPI to specify the number of processes for MPIEXEC; the next option will be the number of processes.")
set(MPIEXEC_POSTFLAGS "-d 1"
  CACHE FILEPATH
  "These flags will come after all flags given to MPIEXEC.")

# Workaround for CMake 3, otherwise FindCUDA sets CUDA_HOST_COMPILER incorrectly.
set(CUDA_HOST_COMPILER "${CMAKE_C_COMPILER}"
  CACHE FILEPATH "Host side compiler used by NVCC")  

# Compile for Kepler compute architecture.
set(CUDA_GPU_ARCH "sm_35" CACHE STRING "" FORCE)

# Change dafault cache entries.
option(DCA_GPU_SUPPORT "Enable GPU support." ON)
option(DCA_PTHREADS    "Enable pthreads"     ON)
