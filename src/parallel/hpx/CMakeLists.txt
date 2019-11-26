# HPX

configure_file("${PROJECT_SOURCE_DIR}/include/dca/config/hpx_defines.hpp.in"
  "${CMAKE_BINARY_DIR}/include/dca/config/hpx_defines.hpp"
  @ONLY)

# working method is

# hpx/CMakeLists.txt for parallel_hpx
#   add_library(parallel_hpx STATIC hpxthread.cpp)
#   hpx_setup_target(parallel_hpx)
#   target_compile_definitions(parallel_hpx INTERFACE DCA_HPX_MAIN)

# dca_testing
#   hpx_setup_target(${name})

add_library(parallel_hpx STATIC hpxthread.cpp)
target_link_libraries(parallel_hpx PRIVATE hpx_wrap)
target_compile_definitions(parallel_hpx INTERFACE DCA_HPX_MAIN)
hpx_setup_target(parallel_hpx)

#target_link_options(parallel_hpx INTERFACE "-Wl,-e,_initialize_main")
