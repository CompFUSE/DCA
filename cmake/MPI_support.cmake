################################################################################
# Checks for MPI C++ compiler wrapper and mpiexec.
#
# TODO: - Make this more robust and more general.
#       - Replace include_directories --> target_include_directories
#       - Replace link_libraries --> target_link_libraries
#
# References: - https://github.com/ALPSCore/ALPSCore
#             - https://github.com/andrealani/ShockFitting
################################################################################

option(DCA_WITH_MPI "Enable MPI support." ON)
set(DCA_HAVE_MPI FALSE CACHE INTERNAL "")

if (DCA_WITH_MPI)
  # Check if CXX compiler already supports MPI.
  include(CheckCXXSourceCompiles)
  check_cxx_source_compiles(
    "#include <mpi.h>
   #include <iostream>
   int main(int argc, char* argv[])
   {
     MPI_Init(&argc, &argv); int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Finalize();
     return 0;
   }"
    CXX_SUPPORTS_MPI)
  
  if (CXX_SUPPORTS_MPI)
    message(STATUS "CXX compiler supports MPI: ${CMAKE_CXX_COMPILER}.")
    set(DCA_HAVE_MPI TRUE)
    
  else()
    # Try to find MPI.
    message(STATUS "CXX compiler does not support MPI. Trying to find MPI.")
    find_package(MPI)
    
    if (${MPI_CXX_FOUND})
      # Check that the versions of compilers are the same.
      execute_process(COMMAND ${MPI_CXX_COMPILER}   "-dumpversion" OUTPUT_VARIABLE mpicxx_version OUTPUT_STRIP_TRAILING_WHITESPACE)
      execute_process(COMMAND ${CMAKE_CXX_COMPILER} "-dumpversion" OUTPUT_VARIABLE cxx_version    OUTPUT_STRIP_TRAILING_WHITESPACE)
      if ("${mpicxx_version}" VERSION_EQUAL "${cxx_version}")
        message(STATUS "Found MPI CXX compiler: ${MPI_CXX_COMPILER}.")
        set(DCA_HAVE_MPI TRUE)
      else()
        message(WARNING "MPI CXX compiler doesn't match CXX compiler.
                       \nMPI support disabled.
                       \nTo enable MPI support set environment variable CXX to mpicxx wrapper
                       \nand run CMake in a fresh build tree.")
        set(DCA_WITH_MPI OFF CACHE BOOL "Enable MPI support." FORCE)
        set(DCA_HAVE_MPI FALSE)
      endif()
      
      list(APPEND CMAKE_CXX_FLAGS ${MPI_CXX_COMPILE_FLAGS})
      include_directories(${MPI_CXX_INCLUDE_PATH} ${MPI_C_INCLUDE_PATH})
      link_libraries(${MPI_CXX_LIBRARIES})
      
    else()
      message(WARNING "MPI not found. MPI support disabled.")
      set(DCA_WITH_MPI OFF CACHE BOOL "Enable MPI support." FORCE)
      set(DCA_HAVE_MPI FALSE)
    endif()
    
  endif()
  
  
  # If MPI is available find MPIEXEC/MPIRUN for execution of tests.
  if (DCA_HAVE_MPI)
    if (NOT ("${MPIEXEC}" STREQUAL "" OR "${MPIEXEC}" STREQUAL "MPIEXEC-NOTFOUND"))
      message(STATUS
        "MPIEXEC already set to ${MPIEXEC}.\n   MPIEXEC_NUMPROC_FLAG is ${MPIEXEC_NUMPROC_FLAG}.")
    else()
      message(STATUS "MPIEXEC not yet set. Trying to find it.")
      find_program(MPIEXEC
        NAMES aprun $ENV{MPIEXEC}
        HINTS ENV MPIEXEC
        )
      if (NOT MPIEXEC_FOUND)
        message(FATAL_ERROR
	  "MPIEXEC not found. Set it manually and reconfigure.")
      else()
        message(STATUS
	  "MPIEXEC set to ${MPIEXEC}.\n   MPIEXEC_NUMPROC_FLAG is ${MPIEXEC_NUMPROC_FLAG}.")
      endif()
    endif()
  endif()
endif()
  
if (DCA_HAVE_MPI)
  dca_add_config_define(DCA_HAVE_MPI)
  set(DCA_PARALLELIZATION_LIBRARY_TYPE "dca::concurrency::MPI_LIBRARY")
  dca_add_config_define(DCA_PARALLELIZATION_LIBRARY_TYPE "dca::concurrency::MPI_LIBRARY")
else()
  dca_add_config_define(DCA_PARALLELIZATION_LIBRARY_TYPE "dca::concurrency::SERIAL_LIBRARY")
  set(DCA_PARALLELIZATION_LIBRARY_TYPE "dca::concurrency::SERIAL_LIBRARY")
endif()
