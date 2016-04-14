#------------------------------------------------------------------------------
# CMake build script for DCA++
#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.0)
project(DCA++ C CXX)
# Disable in-source builds.
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
# include dirs
#------------------------------------------------------------------------------
# For including the generated "dca/config/defines.hpp".
include_directories(${PROJECT_BINARY_DIR}/include)

# For <dca/.../...> #includes.
include_directories(${PROJECT_SOURCE_DIR}/include)

# For "src/.../..." includes.
# Note that these includes should disappear over time since header files should
# be moved into include/dca/.../ and included via #include <dca/.../...>.
include_directories(${PROJECT_SOURCE_DIR})

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
# Modules : for reproducibility
#------------------------------------------------------------------------------
add_subdirectory(modules)
include(modules_init)
get_module_list(MODULE_LIST)
configure_file("${PROJECT_SOURCE_DIR}/modules/modules.cpp.in"
  "${PROJECT_BINARY_DIR}/modules/modules.cpp" @ONLY)

#------------------------------------------------------------------------------
# Git Version : for reproducibility
#------------------------------------------------------------------------------
add_subdirectory(gitVersion)
include(gitVersion_init)
get_git_log(GIT_LOG)  # Executes 'git log -1'.
get_git_status(GIT_STATUS)  # Executes 'git status --porcelain'.
if(GIT_STATUS)
  message(WARNING "Working tree is dirty. Run 'git status' for details.")  
endif()
configure_file("${PROJECT_SOURCE_DIR}/gitVersion/gitVersion.cpp.in"
  "${PROJECT_BINARY_DIR}/gitVersion/gitVersion.cpp" @ONLY)

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
if (DCA_GPU_SUPPORT)
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
  add_subdirectory(src/phys_library/domains/cluster)
  add_subdirectory(src/math_library/geometry_library/vector_operations)
endif()

#------------------------------------------------------------------------------
# Build applications
#------------------------------------------------------------------------------
add_subdirectory(applications)
