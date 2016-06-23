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
set(SPRNG_DIR $ENV{EBROOTSPRNG} 
        CACHE FILEPATH "Path to SPRNG installation directory.")
mark_as_advanced(NFFT_DIR SPGLIB_DIR gtest_DIR SPRNG_DIR)

# The C++ compile wrapper CC already includes and links to these libraries.
# No need to look for them.
set(DCA_HDF5_IMPLICIT   TRUE CACHE INTERNAL "")
set(DCA_LAPACK_IMPLICIT TRUE CACHE INTERNAL "")
set(DCA_FFTW_IMPLICIT   TRUE CACHE INTERNAL "")

# MPIEXEC stuff for executing parallel tests.
# Check whether the system uses the aprun or the srun command. To do this we exploit that executing
# the command that is available results in an error message that contains the command name and
# executing the command that is not available just results in "No such file or directory".
execute_process(COMMAND srun
                RESULT_VARIABLE srun_res
                ERROR_VARIABLE srun_err)
# message ("res = ${srun_res}")
# essage ("err = ${srun_err}")
execute_process(COMMAND aprun
                RESULT_VARIABLE aprun_res
                ERROR_VARIABLE aprun_err)
# message ("res = ${aprun_res}")
# message ("err = ${aprun_err}")

if ("${srun_err}" MATCHES ".*srun.*")
# message ("Use srun.")
set(MPIEXEC "srun"
  CACHE FILEPATH "Executable for running MPI programs.")
elseif ("${aprun_err}" MATCHES ".*aprun.*")
# message ("Use aprun.")
set(MPIEXEC "aprun"
  CACHE FILEPATH "Executable for running MPI programs.")
else ()
message (FATAL_ERROR "Neither aprun nor srun command found.")
endif ()
set(MPIEXEC_NUMPROC_FLAG "-n"
  CACHE FILEPATH "Flag used by MPI to specify the number of processes for
                  MPIEXEC; the next option will be the number of processes.")
mark_as_advanced(MPIEXEC MPIEXEC_NUMPROC_FLAG)
