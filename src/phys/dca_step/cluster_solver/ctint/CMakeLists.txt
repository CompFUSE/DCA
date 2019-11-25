# CT-INT

add_library(ctint
    structs/ctint_matrix_configuration.cpp
    structs/interaction_vertices.cpp
    structs/read_write_configuration.cpp
    structs/solver_configuration.cpp
    )

if (DCA_HAVE_CUDA)
  target_compile_definitions(ctint PUBLIC "-DDCA_HAVE_CUDA")
  target_include_directories(ctint PRIVATE ${MAGMA_INCLUDE_DIR})
  target_link_libraries(ctint PUBLIC ${CUDA_CUDART_LIBRARY};cuda_utils)
endif()


