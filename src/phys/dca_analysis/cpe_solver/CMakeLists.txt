# CPE solver

add_library(cpe_solver STATIC
  basis_functions/delta_function.cpp
  basis_functions/piecewise_linear_function.cpp
  basis_functions/radial_function.cpp
  basis_functions/union_piecewise_linear_and_delta_function.cpp
  basis_functions/union_radial_and_delta_function.cpp)
