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
include(dca_cuda)

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
# Depends on DCA_build_options.
include(dca_compiler_options)

################################################################################
# External libraries
include(dca_external_libs)

################################################################################
# Generate preprocessor definitions header files.
# Depends on DCA_build_options.
dca_write_config_definitions_file()
dca_write_haves_definitions_file()

install(DIRECTORY ${CMAKE_BINARY_DIR}/include/dca/config DESTINATION include/dca)

################################################################################
# Internal libraries
# TODO: When the number of 'internal' libraries grows, we might want to combine these into a few big
#       libraries; e.g. one for math, phys, etc.
add_subdirectory(src/math/geometry/gaussian_quadrature)
add_subdirectory(src/math/random)
add_subdirectory(src/phys_library/DCA_step/cluster_solver/posix_qmci)
add_subdirectory(src/util)

################################################################################
# Common definitions for building DCA++ applications and tests.
# Libraries
set(DCA_LIBRARIES ${DCA_EXTERNAL_LIBS} git_version modules posix_qmci gaussian_quadrature random)

# Includes
set(DCA_INCLUDES ${DCA_EXTERNAL_INCLUDES})

################################################################################
# Testing
option(DCA_WITH_TESTS_FAST "Build DCA++'s fast tests." OFF)
option(DCA_WITH_TESTS_EXTENSIVE "Build DCA++'s extensive tests." OFF)
option(DCA_WITH_TESTS_PERFORMANCE "Build DCA++'s performance tests." OFF)

if (DCA_WITH_TESTS_FAST OR DCA_WITH_TESTS_EXTENSIVE OR DCA_WITH_TESTS_PERFORMANCE)
  include(dca_testing)

  enable_testing()

  add_subdirectory(${gtest_DIR} ${PROJECT_BINARY_DIR}/gtest)

  add_subdirectory(${PROJECT_SOURCE_DIR}/applications/analysis/test)
  add_subdirectory(${PROJECT_SOURCE_DIR}/applications/cluster_solver_check/test)
  add_subdirectory(${PROJECT_SOURCE_DIR}/applications/dca/test)
  add_subdirectory(${PROJECT_SOURCE_DIR}/test)
endif()

################################################################################
# Build applications.
add_subdirectory(applications)
