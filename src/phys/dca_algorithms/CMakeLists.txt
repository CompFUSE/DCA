# DCA(+) algorithms

add_library(dca_algorithms STATIC compute_band_structure.cpp)
dca_gpu_runtime_link(dca_algorithms)
dca_gpu_blas_link(dca_algorithms)
