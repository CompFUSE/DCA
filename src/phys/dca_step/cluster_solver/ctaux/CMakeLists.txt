# CT-AUX

add_library(ctaux STATIC
  domains/feynman_expansion_order_domain.cpp
  domains/hs_field_sign_domain.cpp
  domains/hs_spin_domain.cpp
  domains/hs_vertex_move_domain.cpp
  structs/vertex_singleton.cpp
  walker/ct_aux_walker_tools.cpp)

if (DCA_HAVE_CUDA)
  set_property(TARGET ctaux
               APPEND PROPERTY COMPILE_DEFINITIONS DCA_HAVE_CUDA)
endif ()

if (DCA_HAVE_CUDA)
  CUDA_ADD_LIBRARY(ctaux_walker_kernels
    walker/ct_aux_walker_tools_kernels.cu
    walker/tools/g0_interpolation/g0_interpolation_kernels.cu
    walker/tools/g_matrix_tools/g_matrix_tools_kernels.cu
    walker/tools/n_matrix_tools/n_matrix_tools_kernels.cu)

  target_compile_definitions(ctaux_walker_kernels PRIVATE -DDCA_HAVE_CUDA)
endif()
