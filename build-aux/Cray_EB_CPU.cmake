################################################################################
# Initial cache list for Cray systems where the external libraries are installed
# with EasyBuild.
# CPU version (no GPU support).
# Usage: cmake -C path/to/this/file/Cray_EB_CPU.cmake ...
################################################################################

# Location of external libraries.
# Since they have been installed with EasyBuild, we can use the corresponding
# environment variables set by EB.
set(NFFT_DIR $ENV{EBROOTNFFT}
  CACHE FILEPATH "Path to NFFT installation directory.")
set(SPGLIB_DIR $ENV{EBROOTSPGLIB}
  CACHE FILEPATH "Path to spglib installation directory.")
set(gtest_DIR $ENV{EBROOTGTEST} CACHE FILEPATH "Path to Google Test.")
mark_as_advanced(NFFT_DIR SPGLIB_DIR gtest_DIR)

# The C++ compile wrapper CC already includes and links to these libraries.
# No need to look for them.
set(DCA_HDF5_IMPLICIT   TRUE CACHE INTERNAL "")
set(DCA_LAPACK_IMPLICIT TRUE CACHE INTERNAL "")
set(DCA_FFTW_IMPLICIT   TRUE CACHE INTERNAL "")

# MPIEXEC stuff for executing parallel tests.
set(MPIEXEC "aprun"
  CACHE FILEPATH "Executable for running MPI programs.")
set(MPIEXEC_NUMPROC_FLAG "-n"
  CACHE FILEPATH "Flag used by MPI to specify the number of processes for
                  MPIEXEC; the next option will be the number of processes.")
set(MPIEXEC_POSTFLAGS "-d 1"
  CACHE FILEPATH
  "These flags will come after all flags given to MPIEXEC.")
mark_as_advanced(MPIEXEC MPIEXEC_NUMPROC_FLAG MPIEXEC_POSTFLAGS)
