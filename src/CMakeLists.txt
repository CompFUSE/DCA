add_subdirectory(config)
add_subdirectory(function)
add_subdirectory(io)
add_subdirectory(linalg)
add_subdirectory(math)
add_subdirectory(parallel)
add_subdirectory(phys)
add_subdirectory(profiling)
add_subdirectory(util)

add_library(enumerations STATIC phys/four_point_type.cpp phys/error_computation_type.cpp
            distribution/dist_types.cpp)

