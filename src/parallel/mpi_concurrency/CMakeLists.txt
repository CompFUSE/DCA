# parallel mpi_concurrency
add_library(parallel_mpi_concurrency STATIC mpi_concurrency.cpp mpi_processor_grouping.cpp
            mpi_initializer.cpp mpi_gang.cpp)

if(DCA_HAVE_GPU)
  add_library(kernel_test kernel_test.cu)
  dca_gpu_runtime_link(kernel_test)
  dca_gpu_blas_link(kernel_test)
  dca_gpu_device_link(kernel_test)
  target_link_libraries(kernel_test PUBLIC cuda_utils)
  target_link_libraries(parallel_mpi_concurrency PUBLIC kernel_test)
endif()

# if find_package(MPI) was used we need this, otherwise it does nothing
target_link_libraries(parallel_mpi_concurrency PUBLIC ${MPI_C_LIBRARIES})

