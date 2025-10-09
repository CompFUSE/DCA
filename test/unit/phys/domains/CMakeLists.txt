# domains

add_subdirectory(cluster)
add_subdirectory(time_and_frequency)

add_library(dca_domains INTERFACE)
target_link_libraries(dca_domains INTERFACE  convert_domains
  cluster_domains
  quantum_domains
  time_and_frequency_domains
)
