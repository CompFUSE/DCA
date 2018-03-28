# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)

# Git version
# Note: Pass real source and binary dirs to script.
add_custom_target(git_version_helper
  COMMAND ${CMAKE_COMMAND} -DSCRIPT_SRC_DIR=${PROJECT_SOURCE_DIR}
  -DSCRIPT_BIN_DIR=${PROJECT_BINARY_DIR} -P ${PROJECT_SOURCE_DIR}/cmake/dca_git_version/rebuild.cmake
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

dca_add_library(git_version STATIC ${CMAKE_CURRENT_BINARY_DIR}/git_version.cpp)
add_dependencies(git_version git_version_helper)
if(DCA_HAVE_CUDA)
  target_compile_definitions(git_version PRIVATE -DDCA_HAVE_CUDA)
endif()

# Modules
add_custom_target(modules_helper
  COMMAND ${CMAKE_COMMAND} -DSCRIPT_SRC_DIR=${PROJECT_SOURCE_DIR}
  -DSCRIPT_BIN_DIR=${PROJECT_BINARY_DIR} -P ${PROJECT_SOURCE_DIR}/cmake/dca_modules/rebuild.cmake
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

dca_add_library(modules STATIC ${CMAKE_CURRENT_BINARY_DIR}/modules.cpp)
add_dependencies(modules modules_helper)

# Timer
dca_add_library(timer STATIC timer.cpp)
