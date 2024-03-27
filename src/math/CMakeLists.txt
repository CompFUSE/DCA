add_subdirectory(function_transform)
add_subdirectory(geometry)
add_subdirectory(nfft)
add_subdirectory(random)
add_subdirectory(statistical_testing)

add_library(dca_math INTERFACE)
target_link_libraries(dca_math INTERFACE
  function_transform
  gaussian_quadrature
  tetrahedron_mesh
  nfft
  random
)
