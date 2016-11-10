# DCA(+) algorithms

add_library(dca_algorithms STATIC compute_band_structure.cpp)
target_include_directories(dca_algorithms PUBLIC ${FFTW_INCLUDE_DIR})
