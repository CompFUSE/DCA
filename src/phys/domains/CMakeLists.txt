# Domains

add_subdirectory(cluster)
add_subdirectory(quantum)
add_subdirectory(time_and_frequency)

add_library(convert_domains STATIC convert.cpp)
