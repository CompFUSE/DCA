# Initial cache list for Cray systems where the external libraries have been installed with
# EasyBuild
#
# CPU version (no GPU support)
#
# Usage: cmake -C path/to/this/file/Cray_EB_CPU.cmake ...

# Location of external libraries
# Since they have been installed with EasyBuild, we can use the corresponding environment variables
# that are set by EB.
set(NFFT_DIR $ENV{EBROOTNFFT} CACHE PATH "Path to NFFT installation directory.")
set(SPGLIB_DIR $ENV{EBROOTSPGLIB} CACHE PATH "Path to spglib installation directory.")
set(gtest_DIR $ENV{EBROOTGTEST} CACHE PATH "Path to Google Test.")
set(SPRNG_DIR $ENV{EBROOTSPRNG} CACHE PATH "Path to SPRNG installation directory.")
mark_as_advanced(
  NFFT_DIR
  SPGLIB_DIR
  gtest_DIR
  SPRNG_DIR)

# The C++ compiler wrapper CC already adds the correct include paths for these libraries and
# automatically links to them.
set(DCA_HAVE_FFTW TRUE CACHE INTERNAL "")
set(DCA_HAVE_HDF5 TRUE CACHE INTERNAL "")
set(DCA_HAVE_LAPACK TRUE CACHE INTERNAL "")

# MPIEXEC stuff for executing parallel tests
# If the 'slurm' module is loaded, the command for running MPI programs is 'srun'. Otherwise, check
# whether the 'alps' module is loaded and use 'aprun'.
execute_process(COMMAND modulecmd bash list
  RESULT_VARIABLE res
  ERROR_VARIABLE module_list)

string(FIND ${module_list} "slurm" slurm_found)
if (NOT (${slurm_found} EQUAL -1))
  # Use srun
  set(MPIEXEC "srun"
    CACHE FILEPATH "Executable for running MPI programs.")
  
else()
  # Check for aprun
  string(FIND ${module_list} "alps" alps_found)
  if (NOT (${alps_found} EQUAL -1))
    # Use aprun
    set(MPIEXEC "aprun"
      CACHE FILEPATH "Executable for running MPI programs.")
  else()
    message (FATAL_ERROR "Neither aprun nor srun command found.")
  endif()
  
endif()

set(MPIEXEC_NUMPROC_FLAG "-n"
  CACHE STRING "Flag used by MPI to specify the number of processes for MPIEXEC; the next option will be the number of processes.")
mark_as_advanced(MPIEXEC MPIEXEC_NUMPROC_FLAG)
