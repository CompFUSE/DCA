################################################################################
# Initial cache template for local machines
# Usage: cmake -C path/to/this/file ..
################################################################################

####################################################################
#The user is expected to modify this block
####################################################################
#libraries location.
set(DCA_LIBDIR_NFFT   "" CACHE FILEPATH "NFFT directory") #don't include /lib
set(DCA_LIBDIR_SPGLIB "" CACHE FILEPATH "SPGLIB directory") #don't include /lib
# Compilers for MPI
set(MPIEXEC "" CACHE FILEPATH "Executable for running MPI programs.")
set(CMAKE_CXX_COMPILER CACHE FILEPATH "")           #set to mpi wrapper compiler e.g. mpicxx
set(CMAKE_C_COMPILER CACHE FILEPATH   "")           #set to mpi wrapper compiler e.g. mpicc
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