# parallel mpi_concurrency
add_library(parallel_mpi_concurrency STATIC mpi_concurrency.cpp mpi_processor_grouping.cpp
            mpi_initializer.cpp mpi_gang.cpp)

if(DCA_HAVE_CUDA OR DCA_HAVE_HIP)
  add_library(kernel_test kernel_test.cu)
  if(DCA_HAVE_CUDA)
  set_target_properties(kernel_test PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  set_target_properties(kernel_test PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  target_link_libraries(kernel_test PUBLIC CUDA::cudart gpu_utils)
  target_link_libraries(parallel_mpi_concurrency PUBLIC kernel_test)
  target_compile_definitions(parallel_mpi_concurrency PRIVATE DCA_HAVE_CUDA)
  elseif(DCA_HAVE_HIP)
  set_target_properties(kernel_test PROPERTIES HIP_RESOLVE_DEVICE_SYMBOLS ON)
  set_target_properties(kernel_test PROPERTIES HIP_SEPARABLE_COMPILATION ON)
  target_link_libraries(kernel_test PUBLIC hip::host gpu_utils)
  target_link_libraries(parallel_mpi_concurrency PUBLIC kernel_test)
  target_compile_definitions(parallel_mpi_concurrency PRIVATE DCA_HAVE_HIP)
  endif()
endif()

# if find_package(MPI) was used we need this, otherwise it does nothing
target_link_libraries(parallel_mpi_concurrency PUBLIC ${MPI_C_LIBRARIES})

