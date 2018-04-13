# CT-AUX

add_library(ctint
    domains/ct_int_domains.cpp
    structs/interaction_vertices.cpp
    walker/tools/g0_interpolation.cpp
    walker/tools/d_matrix_builder.cpp
    )

# INTERNAL: consider using --relocatable-device-code true and splitting the .cu files.
if (DCA_HAVE_CUDA)
  CUDA_ADD_LIBRARY(ctint_kernels
      walker/tools/d_matrix_kernels.cu
      walker/gpu_walker_kernels.cu
      OPTIONS  "-DDCA_HAVE_CUDA"
      )
  target_compile_definitions(ctint_kernels PUBLIC  "-DDCA_HAVE_CUDA")
  target_compile_definitions(ctint PRIVATE "-DDCA_HAVE_CUDA")
  target_include_directories(ctint PRIVATE ${MAGMA_INCLUDE_DIR})
  target_link_libraries(ctint PRIVATE ${CUDA_CUDART_LIBRARY})
endif()


