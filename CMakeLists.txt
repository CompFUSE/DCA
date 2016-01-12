################################################################################
# CMake build script for DCA++
################################################################################

cmake_minimum_required(VERSION 3.0)

project(DCA++ C CXX)

# Disable in-source builds.
if (${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
  message(FATAL_ERROR "In-source builds are not permitted.
                      Make a separate folder for building:\n
                      mkdir build; cd build; cmake ..\n
                      Before that, remove the files already created:\n
                      rm -rf CMakeCache.txt CMakeFiles")
endif()


# Include CMake scripts from cmake directory.
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Modules
add_subdirectory(modules)
include(modules_init)
get_module_list(MODULE_LIST)
configure_file("${CMAKE_SOURCE_DIR}/modules/modules.cpp.in"
  "${CMAKE_BINARY_DIR}/modules/modules.cpp" @ONLY)

# Git Version
add_subdirectory(gitVersion)
include(gitVersion_init)
get_git_log(GIT_LOG)  # Executes 'git log -1'.
get_git_status(GIT_STATUS)  # Executes 'git status --porcelain'.
if(GIT_STATUS)
  message(WARNING "Working tree is dirty. Run 'git status' for details.")  
endif()
configure_file("${CMAKE_SOURCE_DIR}/gitVersion/gitVersion.cpp.in"
  "${CMAKE_BINARY_DIR}/gitVersion/gitVersion.cpp" @ONLY)

# DCA++ options
include(DCA_build_options)

# Compiler options and tweaks
include(DCA_compiler_options)

# External libraries
include(DCA_external_libs)

# Common definitions
# NOTE: Needs to be included after DCA_external_libs.
include(DCA_definitions)

# MPI support
option(DCA_MPI_SUPPORT "Enable MPI support." ON)
set(DCA_MPI_AVAILABLE FALSE CACHE INTERNAL "")
if (DCA_MPI_SUPPORT)
  include(MPI_support)
endif()
if (DCA_MPI_AVAILABLE)
  set(DCA_PARALLELIZATION_LIBRARY_TYPE "COMP_LIB::MPI_LIBRARY")
else()
  set(DCA_PARALLELIZATION_LIBRARY_TYPE "COMP_LIB::SERIAL_LIBRARY")
endif()

# GPU support
include(GPU_support)
if (DCA_GPU_SUPPORT)
  add_subdirectory(src/DCA_GPU_routines)
endif()


# Build programs
add_subdirectory(include)

# Build main_dca.cpp (DCA+ loop)
option(DCA_BUILD_DCA "Build main_dca.cpp." ON)
if (DCA_BUILD_DCA)
  add_subdirectory(examples/dca)
endif()

# Build main_analysis.cpp
option(DCA_BUILD_ANALYSIS "Build main_analysis.cpp." ON)
if (DCA_BUILD_ANALYSIS)
  add_subdirectory(examples/analysis)
endif()

# Build check_solver.cpp
option(DCA_BUILD_CHECK_SOLVER "Build check_solver.cpp." ON)
if (DCA_BUILD_CHECK_SOLVER)
  add_subdirectory(examples/check_solver)
endif()


# Testing
option(DCA_TESTS                  "Build DCA++'s tests."                        ON)
option(DCA_TESTS_INCLUDE_ADVANCED "Include time- and resource-consuming tests." ON)
if (DCA_TESTS)
  include(testing)
  add_subdirectory(test)
endif()
