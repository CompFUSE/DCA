################################################################################
# Initial cache list for Cray systems without GPU (Piz Dora)
# Usage: cmake -C path/to/Cray_CPU.cmake ...
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

# Change dafault cache entries.
option(DCA_GPU_SUPPORT "Enable GPU support." OFF)
option(DCA_PTHREADS    "Enable pthreads"     ON)
