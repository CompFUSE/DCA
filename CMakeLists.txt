################################################################################
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#         John Biddiscombe (john.biddiscombe@cscs.ch)
#
# CMake build script for DCA++

cmake_minimum_required(VERSION 3.0)
project(DCA++ C CXX)

################################################################################
# Disable in-source builds.
if (${PROJECT_BINARY_DIR} STREQUAL ${PROJECT_SOURCE_DIR})
  message(FATAL_ERROR "In-source builds are not permitted. Make a separate folder for building:\n
                       mkdir build; cd build; cmake ..\n
                       Before that, remove the files that have already been created:\n
                       rm -rf CMakeCache.txt CMakeFiles")
endif()

################################################################################
# Include CMake scripts from cmake directory.
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

################################################################################
# Load global variables to store preprocessor definitions.
include(dca_defines)

################################################################################
# Check for CXX features that might not be supported.
# If any of these features are not supported by the compiler, then the DCA++ project will not
# compile.
include(dca_cmake_utils)  # for dca_check_cxx_features

dca_check_cxx_features(
  cxx_alias_templates
  cxx_auto_type
  cxx_constexpr
  cxx_decltype
  cxx_decltype_auto
  cxx_decltype_incomplete_return_types
  cxx_default_function_template_args
  cxx_defaulted_functions
  cxx_defaulted_move_initializers
  cxx_generic_lambdas
  cxx_lambdas
  cxx_nullptr
  cxx_range_for
  cxx_rvalue_references
  cxx_static_assert
  cxx_template_template_parameters
  cxx_variadic_templates
)

################################################################################
# Check for libraries and accordingly set DCA_HAVE_XXX variables.

# SPRNG
include(dca_sprng)

# Gnuplot
include(dca_gnuplot)

# MPI
include(dca_mpi)

# Pthreads
include(dca_pthreads)

# CUDA
include(dca_cuda)  # Sets DCA_CUDA_LIBS.

################################################################################
# Include dirs

# For "dca/.../..." #includes.
include_directories(${PROJECT_SOURCE_DIR}/include)

# For src includes, e.g. #include "phys_library/...".
# Note that these includes should disappear over time since header files should be moved into
# include/dca/.../ and included via #include "dca/.../...".
include_directories(${PROJECT_SOURCE_DIR}/src)

# Directory where dca/config/headers.hpp will be created.
include_directories(${CMAKE_BINARY_DIR}/include)

################################################################################
# Modules: for reproducibility
include(dca_modules/init)

get_module_list(MODULE_LIST)

configure_file("${PROJECT_SOURCE_DIR}/src/util/modules.cpp.in"
  "${PROJECT_BINARY_DIR}/src/util/modules.cpp" @ONLY)

################################################################################
# Git Version: for reproducibility
include(dca_git_version/init)

get_git_log(GIT_LOG)        # Executes 'git log -1'.
get_git_status(GIT_STATUS)  # Executes 'git status --porcelain'.

if (GIT_STATUS)
  message(WARNING "Working tree is dirty. Run 'git status' for details.")
endif()

configure_file("${PROJECT_SOURCE_DIR}/src/util/git_version.cpp.in"
  "${PROJECT_BINARY_DIR}/src/util/git_version.cpp" @ONLY)

################################################################################
# Configure DCA++'s applications
include(dca_config)

################################################################################
# Compiler options and tweaks
# Depends on dca_config.
include(dca_compiler_options)

################################################################################
# External libraries
include(dca_external_libs)  # Sets DCA_EXTERNAL_LIBS and DCA_EXTERNAL_INCLUDE_DIRS.

################################################################################
# Generate preprocessor definitions header files.
# Depends on dca_config.
dca_write_config_definitions_file()
dca_write_haves_definitions_file()

install(DIRECTORY ${CMAKE_BINARY_DIR}/include/dca/config DESTINATION include/dca)

################################################################################
# Internal libraries
# TODO: When the number of 'internal' libraries grows, we might want to combine these into a few big
#       libraries; e.g. one for math, phys, etc.
add_subdirectory(src)

if (DCA_HAVE_CUDA)
  add_subdirectory(src/DCA_GPU_routines)
endif()

################################################################################
# Common definitions for building DCA++'s applications and tests.
# Libraries
set(DCA_LIBS
  ${DCA_EXTERNAL_LIBS}
  git_version modules
  posix_qmci
  gaussian_quadrature
  random
  parallel_util
  pthreading
  dca_algorithms
  json_reader json_writer
  hdf5_reader hdf5_writer
  profiling)

if (DCA_HAVE_CUDA)
  list(APPEND DCA_LIBS DCA_GPU_routines ${DCA_CUDA_LIBS})
endif()

# Includes
set(DCA_INCLUDE_DIRS ${DCA_EXTERNAL_INCLUDE_DIRS})

################################################################################
# Testing
option(DCA_WITH_TESTS_FAST "Build DCA++'s fast tests." OFF)
option(DCA_WITH_TESTS_EXTENSIVE "Build DCA++'s extensive tests." OFF)
option(DCA_WITH_TESTS_PERFORMANCE "Build DCA++'s performance tests." OFF)

if (DCA_WITH_TESTS_FAST OR DCA_WITH_TESTS_EXTENSIVE OR DCA_WITH_TESTS_PERFORMANCE)
  include(dca_testing)

  enable_testing()

  if ((CMAKE_BUILD_TYPE STREQUAL "Release") OR
      (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo") OR
      (CMAKE_BUILD_TYPE STREQUAL "MinSizeRel"))
    message(WARNING "Tests are built with NDEBUG defined!")
  endif()

  if (DCA_WITH_TESTS_PERFORMANCE AND NOT (CMAKE_BUILD_TYPE STREQUAL "Release"))
    message(WARNING "Performance tests are only built in Release mode.")
  endif()

  add_subdirectory(${gtest_DIR} ${PROJECT_BINARY_DIR}/gtest)

  add_subdirectory(${PROJECT_SOURCE_DIR}/test)
  add_subdirectory(${PROJECT_SOURCE_DIR}/applications/analysis/test)
  add_subdirectory(${PROJECT_SOURCE_DIR}/applications/cluster_solver_check/test)
  add_subdirectory(${PROJECT_SOURCE_DIR}/applications/dca/test)
endif()

################################################################################
# Build applications.
add_subdirectory(applications)
