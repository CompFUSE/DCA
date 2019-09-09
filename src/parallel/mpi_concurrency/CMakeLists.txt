# parallel mpi_concurrency
add_library(parallel_mpi_concurrency STATIC mpi_concurrency.cpp mpi_processor_grouping.cpp
            mpi_initializer.cpp mpi_gang.cpp)

if(DCA_HAVE_CUDA)
  cuda_add_library(kernel_test kernel_test.cu)
  target_link_libraries(parallel_mpi_concurrency PRIVATE kernel_test)
  target_compile_definitions(parallel_mpi_concurrency PRIVATE DCA_HAVE_CUDA)
endif()

# if find_package(MPI) was used we need this, otherwise it does nothing
target_link_libraries(parallel_mpi_concurrency PUBLIC ${MPI_C_LIBRARIES})

