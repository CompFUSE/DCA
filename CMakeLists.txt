################################################################################
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#         John Biddiscombe (john.biddiscombe@cscs.ch)
#
# CMake build script for DCA++

cmake_minimum_required(VERSION 3.3)
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
  cxx_variadic_templates)

################################################################################
# Modules: for reproducibility
include(dca_modules/init)
get_module_list(MODULE_LIST)

configure_file("${PROJECT_SOURCE_DIR}/src/util/modules.cpp.in"
  "${PROJECT_BINARY_DIR}/src/util/modules.cpp" @ONLY)

################################################################################
# Git version: for reproducibility
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
# Generate preprocessor definitions header files.
# Depends on dca_config.
dca_write_config_definitions_file()
dca_write_haves_definitions_file()

install(DIRECTORY ${CMAKE_BINARY_DIR}/include/dca/config DESTINATION include/dca)

################################################################################
# External libraries
include(dca_external_libs)  # Sets DCA_EXTERNAL_LIBS and DCA_EXTERNAL_INCLUDE_DIRS.

################################################################################
# Include directories
# For "dca/.../..." #includes.
include_directories(${PROJECT_SOURCE_DIR}/include)
# Directory where dca/config/headers.hpp will be created.
include_directories(${CMAKE_BINARY_DIR}/include)

set(DCA_INCLUDE_DIRS
  ${DCA_EXTERNAL_INCLUDE_DIRS})

################################################################################
# Libraries
add_subdirectory(src)

set(DCA_LIBS
  git_version
  modules
  ctaux
  posix_qmci
  ss_ct_hyb
  function_transform
  gaussian_quadrature
  nfft
  tetrahedron_mesh
  ${DCA_RNG_LIBRARY}
  parallel_util
  pthreading
  dca_algorithms
  json
  hdf5
  profiling
  function
  convert_domains
  cluster_domains
  quantum_domains
  time_and_frequency_domains
  symmetrization
  coarsegraining
  ${DCA_EXTERNAL_LIBS})

if (DCA_HAVE_CUDA)
  list(APPEND DCA_LIBS
    ctaux_walker_kernels
    lapack_kernels
    blas_kernels
    cuda_utils
    ${DCA_CUDA_LIBS})
endif()

################################################################################
# Testing
option(DCA_WITH_TESTS_FAST "Build DCA++'s fast tests." OFF)
option(DCA_WITH_TESTS_EXTENSIVE "Build DCA++'s extensive tests." OFF)
option(DCA_WITH_TESTS_PERFORMANCE "Build DCA++'s performance tests. (only in release mode)" OFF)

set(TEST_RUNNER "" CACHE STRING "Command for executing (MPI) programs.")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "Flag used by TEST_RUNNER to specify the number of processes.")

set(gtest_dir "${PROJECT_SOURCE_DIR}/libs/googletest-1.8.0" CACHE PATH "Path to Google Test.")
mark_as_advanced(gtest_dir)

if (DCA_WITH_TESTS_FAST OR DCA_WITH_TESTS_EXTENSIVE OR DCA_WITH_TESTS_PERFORMANCE)
  include(dca_testing)
  enable_testing()

  if ((CMAKE_BUILD_TYPE STREQUAL "Release") OR
      (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo") OR
      (CMAKE_BUILD_TYPE STREQUAL "MinSizeRel"))
    message(WARNING "Tests are built with NDEBUG defined!")
  endif()

  if (DCA_WITH_TESTS_PERFORMANCE AND NOT (CMAKE_BUILD_TYPE STREQUAL "Release"))
    message(WARNING "Performance tests are only built in release mode.")
  endif()

  if (DCA_HAVE_MPI AND (TEST_RUNNER STREQUAL ""))
    message(FATAL_ERROR "TEST_RUNNER needs to be set to the command for executing MPI programs, e.g. mpiexec, mpirun, aprun or srun.")
  endif()
  
  add_subdirectory(${gtest_dir})
  add_subdirectory(${PROJECT_SOURCE_DIR}/test)
  add_subdirectory(${PROJECT_SOURCE_DIR}/applications/cluster_solver_check/test)
endif()

################################################################################
# Build applications.
add_subdirectory(applications)
