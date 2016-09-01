################################################################################
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# This cmake script is executed every time 'make' is executed. It checks the list of currently
# loaded modules (module list). If the list changes, modules.cpp will be reconfigured, which
# triggers the recompilation of libmodules and 'make clean' is executed, which triggers the
# recompilation of all targets.
#
# The param SCRIPT_BIN_DIR is passed with the correct binary path because when this script runs,
# PROJECT_BINARY_DIR is not set as expected.

execute_process(
  COMMAND modulecmd bash list
  WORKING_DIRECTORY "${SCRIPT_SRC_DIR}"
  RESULT_VARIABLE res
  ERROR_VARIABLE err)

set(MODULE_LIST "${err}")
string(REPLACE "\n" "\\n" MODULE_LIST "${MODULE_LIST}")

if (EXISTS "${SCRIPT_BIN_DIR}/src/util/module_list.txt")
  file(READ "${SCRIPT_BIN_DIR}/src/util/module_list.txt" module_list_txt)

  if (NOT (err STREQUAL module_list_txt))
    file(WRITE "${SCRIPT_BIN_DIR}/src/util/module_list.txt" "${err}")
    configure_file("${SCRIPT_SRC_DIR}/src/util/modules.cpp.in" "${SCRIPT_BIN_DIR}/src/util/modules.cpp" @ONLY)

    # Execute 'make clean'.
    execute_process(
      COMMAND make clean
      WORKING_DIRECTORY "${SCRIPT_BIN_DIR}")  # .../build
  endif()

else()
  file(WRITE "${SCRIPT_BIN_DIR}/modules/module_list.txt" "${err}")
  configure_file("${SCRIPT_SRC_DIR}/src/util/modules.cpp.in" "${SCRIPT_BIN_DIR}/src/util/modules.cpp" @ONLY)

  execute_process(
    COMMAND make clean
    WORKING_DIRECTORY "${SCRIPT_BIN_DIR}")
endif()
