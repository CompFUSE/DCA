# CT-AUX

if (DCA_HAVE_MULDEFS)
add_library(ctaux
  STATIC
  domains/feynman_expansion_order_domain.cpp
  domains/hs_field_sign_domain.cpp
  domains/hs_spin_domain.cpp
  domains/hs_vertex_move_domain.cpp
  structs/vertex_singleton.cpp
  walker/ct_aux_walker_tools.cpp)

else()
  add_library(ctaux
    SHARED
    domains/feynman_expansion_order_domain.cpp
    domains/hs_field_sign_domain.cpp
    domains/hs_spin_domain.cpp
    domains/hs_vertex_move_domain.cpp
    structs/vertex_singleton.cpp
    walker/ct_aux_walker_tools.cpp)

  target_link_libraries(ctaux PUBLIC ${LAPACK_LIBRARIES})
endif()
