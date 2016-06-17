################################################################################
# applications/CMakeLists.txt
################################################################################

# Directory of tight_binding_on_2D_square_lattice.hpp etc.
set(DCA_LATTICE_MODELS_FOLDER "${CMAKE_CURRENT_SOURCE_DIR}/lattice_models")

# Build main_dca.cpp (DCA+ loop)
option(DCA_BUILD_DCA "Build main_dca.cpp." ON)
add_subdirectory(dca)

# Build main_analysis.cpp
option(DCA_BUILD_ANALYSIS "Build main_analysis.cpp." ON)
add_subdirectory(analysis)

# Build cluster_solver_check.cpp
option(DCA_BUILD_CLUSTER_SOLVER_CHECK "Build cluster_solver_check." OFF)
add_subdirectory(cluster_solver_check)
