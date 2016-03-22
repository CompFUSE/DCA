# File: cmake/modules_init.cmake
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# Defines a function to get the list of loaded modules (module list).

function(get_module_list _module_list)
  
  execute_process(
    COMMAND
    modulecmd
    bash
    list
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    RESULT_VARIABLE res
    ERROR_VARIABLE err
    # ERROR_STRIP_TRAILING_WHITESPACE
    )

  file(WRITE "${PROJECT_BINARY_DIR}/modules/module_list.txt" "${err}")
  string(REPLACE "\n" "\\n" err "${err}")
  set(${_module_list} "${err}" PARENT_SCOPE)

endfunction()
