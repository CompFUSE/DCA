# File: cmake/gitVersion_rebuild.cmake
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# This cmake script is executed every time 'make' is executed.
# It checks the list of loaded modules (module list).
# If the list changes, modules.cpp will be reconfigured which triggers the
# recompilation of libversion and 'make clean' is executed which triggers
# the recompilation of all targets.

execute_process(
  COMMAND
  modulecmd
  bash
  list
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"  # .../build
  RESULT_VARIABLE   res
  ERROR_VARIABLE    err
  # ERROR_STRIP_TRAILING_WHITESPACE
  )

set(MODULE_LIST "${err}")
string(REPLACE "\n" "\\n" MODULE_LIST "${MODULE_LIST}")

if (EXISTS "${CMAKE_BINARY_DIR}/modules/module_list.txt")
  file(READ "${CMAKE_BINARY_DIR}/modules/module_list.txt" module_list_txt)
  
  if (NOT (err STREQUAL module_list_txt))
    file(WRITE "${CMAKE_BINARY_DIR}/modules/module_list.txt" "${err}")
    configure_file("${CMAKE_BINARY_DIR}/../modules/modules.cpp.in" "${CMAKE_BINARY_DIR}/modules/modules.cpp" @ONLY)
  
    # Execute 'make clean'.
    execute_process(
      COMMAND
      make
      clean
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
      )
  endif()
  
else()
  file(WRITE "${CMAKE_BINARY_DIR}/modules/module_list.txt" "${err}")
  configure_file("${CMAKE_BINARY_DIR}/../modules/modules.cpp.in" "${CMAKE_BINARY_DIR}/modules/modules.cpp" @ONLY)
  
  execute_process(
    COMMAND
    make
    clean
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    )
  
endif()
