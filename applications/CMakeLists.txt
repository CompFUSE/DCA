################################################################################
# applications/CMakeLists.txt
################################################################################

# Copy type_definitions.h only during the initial configuration, otherwise
# changes to the 'build copy' might be overwritten when reconfiguring.
<<<<<<< 03b20f4141d8d2ccba1de25a7f16f0966a44c6ee

if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/model_type.hpp")
  file(COPY model_type.hpp DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/lattice_types.hpp")
    file(COPY lattice_types.hpp DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()
#TODO add rng definition

set(DCA_APPLICATION_FOLDER ${CMAKE_CURRENT_BINARY_DIR}) #dir of model_type.hpp

# Build main_dca.cpp (DCA+ loop)
option(DCA_BUILD_DCA "Build main_dca.cpp." ON)
add_subdirectory(dca)

# Build main_analysis.cpp
option(DCA_BUILD_ANALYSIS "Build main_analysis.cpp." ON)
add_subdirectory(analysis)

# Build cluster_solver_vs_ED.cpp
option(DCA_BUILD_CLUSTER_SOLVER_CHECK "Build cluster_solver_check." OFF)
add_subdirectory(cluster_solver_check)
