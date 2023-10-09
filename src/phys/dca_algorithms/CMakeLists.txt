# DCA(+) algorithms

# add_library(dca_algorithms INTERFACE)

# target_sources(dca_algorithms
#   include/dca/phys/dca_algorithms
  
add_library(dca_algorithms STATIC compute_band_structure.cpp)
dca_gpu_runtime_link(dca_algorithms)
