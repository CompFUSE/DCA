# CT-AUX

add_library(ctaux STATIC
  domains/feynman_expansion_order_domain.cpp
  domains/hs_field_sign_domain.cpp
  domains/hs_spin_domain.cpp
  domains/hs_vertex_move_domain.cpp
  structs/vertex_singleton.cpp
  walker/ct_aux_walker_tools.cpp)

dca_gpu_runtime_link(ctaux)
