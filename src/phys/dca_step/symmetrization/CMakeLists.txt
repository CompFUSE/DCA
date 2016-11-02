# Symmetrization

add_library(symmetrization STATIC symmetrize_single_particle_function.cpp)
include_directories(symmetrization PUBLIC ${FFTW_INCLUDE_DIR})
