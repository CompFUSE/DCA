# parallel

if (DCA_HAVE_MPI)
  add_subdirectory(mpi_concurrency)
endif(DCA_HAVE_MPI)

add_subdirectory(no_concurrency)
add_subdirectory(no_threading)
if (DCA_HAVE_HPX)
  add_subdirectory(hpx)
else()
  add_subdirectory(stdthread)
endif()
add_subdirectory(util)
