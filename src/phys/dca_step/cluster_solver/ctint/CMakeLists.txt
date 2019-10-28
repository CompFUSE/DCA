# CT-INT

if (DCA_HAVE_CUDA)
  set(CTINT_CUDA_CPP_FILES walker/tools/d_matrix_builder_gpu.cpp)
endif ()

add_library(ctint
    domains/ct_int_domains.cpp
    structs/interaction_vertices.cpp
    structs/read_write_configuration.cpp
    walker/tools/g0_interpolation.cpp
    walker/tools/d_matrix_builder.cpp
    ${CTINT_CUDA_CPP_FILES}
    )

if (DCA_HAVE_CUDA)
  target_compile_definitions(ctint PUBLIC "-DDCA_HAVE_CUDA")
  target_include_directories(ctint PRIVATE ${MAGMA_INCLUDE_DIR})
  target_link_libraries(ctint PUBLIC ${CUDA_CUDART_LIBRARY};cuda_utils)
endif()


