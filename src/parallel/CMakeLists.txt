# parallel

if (DCA_THREADING_LIBRARY STREQUAL POSIX)
  add_subdirectory(pthreading)
endif()

if (DCA_THREADING_LIBRARY STREQUAL STDTHREAD)
  add_subdirectory(stdthread)
endif()

if (DCA_HAVE_MPI)
  add_subdirectory(mpi_concurrency)
endif(DCA_HAVE_MPI)

add_subdirectory(no_concurrency)
add_subdirectory(no_threading)
add_subdirectory(util)
