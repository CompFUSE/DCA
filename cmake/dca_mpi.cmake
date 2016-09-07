################################################################################
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# Checks for MPI and accordingly sets DCA_HAVE_MPI.

set(DCA_HAVE_MPI FALSE CACHE INTERNAL "")

# Check if CXX compiler supports MPI.
include(CheckCXXSourceCompiles)

check_cxx_source_compiles(
  "#include <mpi.h>
   int main(int argc, char** argv)
   {
     MPI_Init(&argc, &argv);
     MPI_Finalize();
     return 0;
   }"
  CXX_SUPPORTS_MPI)

if (CXX_SUPPORTS_MPI)
  set(DCA_HAVE_MPI TRUE CACHE INTERNAL "")
  dca_add_haves_define(DCA_HAVE_MPI)

  # Check if MPIEXEC and MPI_NUMPROC_FLAG are set.
  # They are needed to run some of the tests.
  if (NOT MPIEXEC OR NOT MPIEXEC_NUMPROC_FLAG)
    message(FATAL_ERROR "MPIEXEC and MPIEXEC_NUMPROC_FLAG must be set.")
  endif()
endif()
