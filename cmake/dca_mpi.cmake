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
else()
  # if MPICC is not the default compiler, try finding MPI
  # using the usual CMake find_package mechanism
  find_package(MPI QUIET)
  if (MPI_FOUND)
    set(DCA_HAVE_MPI TRUE CACHE INTERNAL "")
    dca_add_haves_define(DCA_HAVE_MPI)
    include_directories(${MPI_C_INCLUDE_PATH})
  endif()
endif()
