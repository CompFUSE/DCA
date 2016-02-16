################################################################################
# CMake build script for DCA++
################################################################################

cmake_minimum_required(VERSION 3.0)
project(DCA++ C CXX)

MESSAGE("project dir ${CMAKE_SOURCE_DIR}")


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
include(MPI_support)

# GPU support
include(GPU_support)
if (DCA_GPU_SUPPORT)
  add_subdirectory(src/DCA_GPU_routines)
endif()

# Testing
include(testing)

# Build applications
add_subdirectory(applications)
