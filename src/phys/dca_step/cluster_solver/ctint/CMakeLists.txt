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

# INTERNAL: consider using --relocatable-device-code true and splitting the .cu files.
if (DCA_HAVE_CUDA)
  CUDA_ADD_LIBRARY(ctint_kernels
      walker/walker_kernels.cu
      walker/tools/d_matrix_kernels.cu
      )
  target_compile_definitions(ctint_kernels PUBLIC  "-DDCA_HAVE_CUDA")

  target_compile_definitions(ctint PRIVATE "-DDCA_HAVE_CUDA")
  target_include_directories(ctint PRIVATE ${MAGMA_INCLUDE_DIR})
  target_link_libraries(ctint PRIVATE ${CUDA_CUDART_LIBRARY})
endif()


