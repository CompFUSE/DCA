# CT-INT

if (DCA_HAVE_CUDA)
  set(CTINT_CUDA_CPP_FILES walker/tools/d_matrix_builder_gpu.cpp)
endif ()

add_library(ctint
    domains/ct_int_domains.cpp
    structs/interaction_vertices.cpp
    walker/tools/g0_interpolation.cpp
    walker/tools/d_matrix_builder.cpp
    structs/read_write_configuration.cpp
    ${CTINT_CUDA_CPP_FILES}
    )

if (DCA_HAVE_CUDA)
  CUDA_ADD_LIBRARY(ctint_kernels
      walker/walker_kernels.cu
      walker/tools/d_matrix_kernels.cu
      device_helper/ctint_helper.cu
      walker/tools/g0_interpolation_gpu.cu
      ../shared_tools/cluster_helper.cu
      )

  set_target_properties(ctint_kernels PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  target_link_libraries(ctint_kernels mc_tools_kernels)
  target_compile_definitions(ctint_kernels PUBLIC  "-DDCA_HAVE_CUDA")



  target_compile_definitions(ctint PUBLIC "-DDCA_HAVE_CUDA")
  target_include_directories(ctint PRIVATE ${MAGMA_INCLUDE_DIR})
  target_link_libraries(ctint PUBLIC ${CUDA_CUDART_LIBRARY};cuda_utils)
endif()


