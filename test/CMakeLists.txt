# test

# Unit tests
add_subdirectory(unit)

# Performance tests
# Only build the performance tests in Release mode.
if (CMAKE_BUILD_TYPE STREQUAL "Release")
  add_subdirectory(performance)
endif()
