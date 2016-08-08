# test

# Only build performance tests in Release mode.
if (CMAKE_BUILD_TYPE STREQUAL "Release")
  add_subdirectory(performance)
endif()

add_subdirectory(unit)
