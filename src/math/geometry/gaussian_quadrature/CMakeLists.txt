# Gaussian quadrature

add_library(gaussian_quadrature STATIC compute_weights_and_abscissas.cpp)
target_include_directories(gaussian_quadrature PUBLIC ${SIMPLEX_GM_RULE_INCLUDE_DIR})
target_link_libraries(gaussian_quadrature PRIVATE ${SIMPLEX_GM_RULE_LIBRARY})
