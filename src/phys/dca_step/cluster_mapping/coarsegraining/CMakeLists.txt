# Coarse-graining

add_library(coarsegraining
  STATIC
  coarsegrain_domain_names.cpp
  tetrahedron_routines_harmonic_function.cpp)
target_include_directories(coarsegraining PUBLIC ${FFTW_INCLUDE_DIR})
