# CT-INT

if (DCA_HAVE_GPU)
  set(CTINT_CUDA_CPP_FILES walker/tools/d_matrix_builder_gpu.cpp)
endif ()

add_library(ctint
    structs/interaction_vertices.cpp
    structs/ctint_matrix_configuration.cpp
    structs/read_write_configuration.cpp
    structs/solver_configuration.cpp
    walker/tools/d_matrix_builder.cpp
    structs/utils.cpp
    ${CTINT_CUDA_CPP_FILES}
    )

target_link_libraries(ctint PUBLIC time_and_frequency_domains dca_io json)

if (DCA_HAVE_GPU)
  target_include_directories(ctint PRIVATE ${MAGMA_INCLUDE_DIR})  
  target_link_libraries(ctint PUBLIC gpu_utils)
  dca_gpu_runtime_link(ctint)
endif()
