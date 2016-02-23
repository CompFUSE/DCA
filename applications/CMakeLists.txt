################################################################################
# applications/CMakeLists.txt
################################################################################

# Copy type_definitions.h only during the initial configuration, otherwise
# changes to the 'build copy' might be overwritten when reconfiguring.
configure_file(type_definitions.h "${CMAKE_BINARY_DIR}/applications/type_definitions.h" @ONLY)
#if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/type_definitions.h")
#  file(COPY type_definitions.h DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
#endif()

# Build main_dca.cpp (DCA+ loop)
option(DCA_BUILD_DCA "Build main_dca.cpp." ON)
add_subdirectory(dca)

# Build main_analysis.cpp
option(DCA_BUILD_ANALYSIS "Build main_analysis.cpp." ON)
add_subdirectory(analysis)

# Build cluster_solver_vs_ED.cpp
option(DCA_BUILD_CLUSTER_SOLVER_CHECK "Build cluster_solver_check." OFF)
add_subdirectory(cluster_solver_check)
