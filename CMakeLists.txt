#------------------------------------------------------------------------------
# CMake build script for DCA++
#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.0)
project(DCA++ C CXX)

#------------------------------------------------------------------------------
# Disable in-source builds.
#------------------------------------------------------------------------------
if (${PROJECT_BINARY_DIR} STREQUAL ${PROJECT_SOURCE_DIR})
  message(FATAL_ERROR "In-source builds are not permitted.
                      Make a separate folder for building:\n
                      mkdir build; cd build; cmake ..\n
                      Before that, remove the files already created:\n
                      rm -rf CMakeCache.txt CMakeFiles")
endif()

#------------------------------------------------------------------------------
# Include CMake scripts from cmake directory.
#------------------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#------------------------------------------------------------------------------
# load project macros
#------------------------------------------------------------------------------
include(DCA_macros)

#------------------------------------------------------------------------------
# add the project source dir to a global #define : we use this as the
# root path for getting files for/from testing
#------------------------------------------------------------------------------
dca_add_config_define(DCA_SOURCE_DIRECTORY \"${PROJECT_SOURCE_DIR}\")

#------------------------------------------------------------------------------
# Include dirs
#------------------------------------------------------------------------------
# For "dca/.../..." #includes.
include_directories(${PROJECT_SOURCE_DIR}/include)

# For src includes, e.g. #include "phys_library/...".
# Note that these includes should disappear over time since header files should
# be moved into include/dca/.../ and included via #include "dca/.../...".
include_directories(${PROJECT_SOURCE_DIR}/src)

# Directory where dca/config/defines.hpp will be created.
include_directories(${CMAKE_BINARY_DIR})

#------------------------------------------------------------------------------
# Check for CXX features that might not be supported 
# if any of these features are not supported by the compiler, then
# the DCA+ project will not compile, so check them
#------------------------------------------------------------------------------
include(DCA_cmake_features)

DCA_check_cxx_features(
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

#------------------------------------------------------------------------------
# Include CMake scripts from cmake directory.
#------------------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#------------------------------------------------------------------------------
# Modules : for reproducibility
#------------------------------------------------------------------------------

include(modules_init)
get_module_list(MODULE_LIST)
configure_file("${PROJECT_SOURCE_DIR}/src/util/modules.cpp.in"
  "${PROJECT_BINARY_DIR}/src/util/modules.cpp" @ONLY)

#------------------------------------------------------------------------------
# Git Version : for reproducibility
#------------------------------------------------------------------------------
include(git_version_init)
get_git_log(GIT_LOG)  # Executes 'git log -1'.
get_git_status(GIT_STATUS)  # Executes 'git status --porcelain'.
if(GIT_STATUS)
  message(WARNING "Working tree is dirty. Run 'git status' for details.")  
endif()
configure_file("${PROJECT_SOURCE_DIR}/src/util/git_version.cpp.in"
  "${PROJECT_BINARY_DIR}/src/util/git_version.cpp" @ONLY)

#------------------------------------------------------------------------------
# DCA++ options
#------------------------------------------------------------------------------
include(DCA_build_options)

#------------------------------------------------------------------------------
# Compiler options and tweaks
#------------------------------------------------------------------------------
include(DCA_compiler_options)

#------------------------------------------------------------------------------
# External libraries
#------------------------------------------------------------------------------
include(DCA_external_libs)

#------------------------------------------------------------------------------
# Internal libraries
# TODO: When the number of 'internal' libraries grows, we might want to combine
#       these into a few big libraries; e.g. one for comp, math and phys,
#       respectively.
#------------------------------------------------------------------------------
add_subdirectory(src/math/geometry/gaussian_quadrature)
add_subdirectory(src/util)
add_subdirectory(src/phys_library/DCA_step/cluster_solver/posix_qmci)

#------------------------------------------------------------------------------
# Common definitions
# NOTE: Needs to be included after DCA_external_libs.
#------------------------------------------------------------------------------
include(DCA_definitions)

#------------------------------------------------------------------------------
# MPI support
#------------------------------------------------------------------------------
include(MPI_support)

#------------------------------------------------------------------------------
# GPU support
#------------------------------------------------------------------------------
include(GPU_support)
if (DCA_WITH_CUDA)
  add_subdirectory(src/DCA_GPU_routines)
endif()

#------------------------------------------------------------------------------
# Global defines.hpp file
#------------------------------------------------------------------------------
write_global_definitions_file("include/dca/config/defines.hpp")

#------------------------------------------------------------------------------
# Testing
#------------------------------------------------------------------------------
include(testing)

if (DCA_TESTS)
  add_subdirectory(src/comp_library/function_library/test)
  add_subdirectory(src/math_library/geometry_library/vector_operations)
endif()

#------------------------------------------------------------------------------
# Configure the header containing all compile definitions
#------------------------------------------------------------------------------
get_property(DCA_CONFIG_DEFINITIONS_VAR GLOBAL PROPERTY DCA_CONFIG_DEFINITIONS)

list(SORT DCA_CONFIG_DEFINITIONS_VAR)
list(REMOVE_DUPLICATES DCA_CONFIG_DEFINITIONS_VAR)
list(REMOVE_ITEM DCA_CONFIG_DEFINITIONS_VAR "")

set(dca_config_defines "// Generated from CMake definitons\n\n")
foreach(def ${DCA_CONFIG_DEFINITIONS_VAR})
  set(dca_config_defines "${dca_config_defines}#define ${def} ${${def}_define}\n")
endforeach()

# Generate a defines.hpp to be used in the build directory ...
set(DCA_DEFINES_PREFIX ${DCA_BUILD_PREFIX})
configure_file("${PROJECT_SOURCE_DIR}/cmake/templates/config_defines.hpp.in"
               "${CMAKE_BINARY_DIR}/dca/config/defines.hpp"
               @ONLY)

# Generate a defines.hpp to be used in the install directory ...
set(DCA_DEFINES_PREFIX ${DCA_PREFIX})
configure_file("${PROJECT_SOURCE_DIR}/cmake/templates/config_defines.hpp.in"
               "${CMAKE_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/dca/config/defines.hpp"
               @ONLY)

#------------------------------------------------------------------------------
# Build applications
#------------------------------------------------------------------------------
add_subdirectory(applications)
