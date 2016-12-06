# Initial cache list for Cray systems where the external libraries have been installed with
# EasyBuild
#
# Multicore version (no GPU support)
#
# Usage: cmake -C /path/to/this/file [other options] /path/to/source/dir

# Location of SPRNG
set(SPRNG_DIR $ENV{EBROOTSPRNG} CACHE PATH "Path to SPRNG installation directory.")

# The C++ compiler wrapper CC already adds the correct include paths for these libraries and
# automatically links to them.
set(DCA_HAVE_FFTW TRUE CACHE INTERNAL "")
set(DCA_HAVE_HDF5 TRUE CACHE INTERNAL "")
set(DCA_HAVE_LAPACK TRUE CACHE INTERNAL "")

# Find the command to run an application, which we will call TEST_RUNNER.
# If the 'slurm' module is loaded, the command is 'srun'. Otherwise, check whether the 'alps' module
# is loaded and use 'aprun'.
execute_process(COMMAND modulecmd bash list
  RESULT_VARIABLE res
  ERROR_VARIABLE module_list)

string(FIND ${module_list} "slurm" slurm_found)
if (NOT (${slurm_found} EQUAL -1))
  # Use srun
  set(TEST_RUNNER "srun"
    CACHE FILEPATH "Command to run an application.")

else()
  # Check for aprun
  string(FIND ${module_list} "alps" alps_found)
  if (NOT (${alps_found} EQUAL -1))
    # Use aprun
    set(TEST_RUNNER "aprun"
      CACHE FILEPATH "Command to run an application.")
  else()
    message (FATAL_ERROR "Neither aprun nor srun command found.")
  endif()

endif()

set(MPIEXEC_NUMPROC_FLAG "-n"
  CACHE STRING "Flag used by TEST_RUNNER to specify the number of processes.")
