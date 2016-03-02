
################################################################################
# Initial cache template for local machines
# Usage: 
# Set CXX and CC to the desired compiler in the environment
# cmake -C path/to/this/file ..
################################################################################

####################################################################
#The user is expected to modify this block
####################################################################
#libraries location.
set(DCA_LIBDIR_NFFT   "/path/to/nfft_root" CACHE FILEPATH "NFFT directory") #don't include /lib
set(DCA_LIBDIR_SPGLIB "/path_to/spglib_root" CACHE FILEPATH "SPGLIB directory") #don't include /lib
# Compilers for MPI
set(MPIEXEC "/path/to/mpirun" CACHE FILEPATH "Executable for running MPI programs.")
####################################################################

#google tests directory
set(GTEST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../libs/gmock-1.7.0/gtest/" CACHE FILEPATH "folder with gtest")
#options and flags
option(DCA_GPU_SUPPORT "Disable GPU support." OFF)
# MPIEXEC stuff for executing parallel tests.
# Note: CXX and CC must be set as environment variables to the corresponding MPI
#       C++ and C compiler wrappers.
set(MPIEXEC_NUMPROC_FLAG "-np" CACHE FILEPATH "Flag used by MPI to specify the number of processes for MPIEXEC; the next option will be the number of processes.")

mark_as_advanced(GTEST_DIR MPIEXEC_NUMPROC_FLAG)