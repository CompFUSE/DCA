################################################################################
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# Defines a function to get a list of the currently loaded modules (module list).

function(get_module_list _module_list)
    execute_process(
    COMMAND modulecmd bash list
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    RESULT_VARIABLE res
    ERROR_VARIABLE err)

  file(WRITE "${PROJECT_BINARY_DIR}/src/util/module_list.txt" "${err}")
  string(REPLACE "\n" "\\n" err "${err}")
  set(${_module_list} "${err}" PARENT_SCOPE)
endfunction()
