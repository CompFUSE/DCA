# File: cmake/gitVersion_init.cmake
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# Defines functions to get the git log (git log -1) and the git status
# (git status --porcelain).

# Check git log
function(get_git_log _git_log)

  execute_process(
    COMMAND
    git --git-dir ${PROJECT_SOURCE_DIR}/.git --work-tree ${PROJECT_SOURCE_DIR}
    log -1
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    RESULT_VARIABLE res
    OUTPUT_VARIABLE out
    )

  file(WRITE "${CMAKE_BINARY_DIR}/gitVersion/git_log.txt" "${out}")
  string(REPLACE "\n" "\\n" out "${out}")
  set(${_git_log} "${out}" PARENT_SCOPE)

endfunction()

# Check git status
function(get_git_status _git_status)

  execute_process(
    COMMAND
    git --git-dir ${PROJECT_SOURCE_DIR}/.git --work-tree ${PROJECT_SOURCE_DIR}
    status --porcelain
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    RESULT_VARIABLE res
    OUTPUT_VARIABLE out
    )

  file(WRITE "${CMAKE_BINARY_DIR}/gitVersion/git_status.txt" "${out}")
  string(REPLACE "\n" "\\n" out "${out}")
  set(${_git_status} "${out}" PARENT_SCOPE)

endfunction()
