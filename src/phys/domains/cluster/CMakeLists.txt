# Cluster domains

add_library(cluster_domains STATIC
  cluster_definitions.cpp
  momentum_exchange_domain.cpp
  symmetries/symmetry_operations/group_action.cpp)

dca_gpu_runtime_link(cluster_domains)
