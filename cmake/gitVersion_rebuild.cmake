# File: cmake/gitVersion_rebuild.cmake
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# This cmake script is executed every time 'make' is executed.
# It checks the last git log (git log -1) and git status
# (git status --porcelain applications src testing).
# If either of them changes, gitVersion.cpp will be reconfigured which triggers the
# recompilation of libgitVersion.

# SCRIPT_SRC_DIR is passed in from the custom rule
# as the PROJECT_SOURCE_DIR is not valid when run as a rule

# Check git log
execute_process(
  COMMAND
  git --git-dir ${SCRIPT_SRC_DIR}/.git --work-tree ${SCRIPT_SRC_DIR}
  log -1
  WORKING_DIRECTORY "${SCRIPT_BIN_DIR}"  # .../build
  RESULT_VARIABLE res
  OUTPUT_VARIABLE out
  )

set(GIT_LOG "${out}")
# Replace newlines with characters "\n" for use in std::string.
string(REPLACE "\n" "\\n" GIT_LOG "${GIT_LOG}")

set(LOG_CHANGED FALSE)

if (EXISTS "${SCRIPT_BIN_DIR}/gitVersion/git_log.txt")
  file(READ "${SCRIPT_BIN_DIR}/gitVersion/git_log.txt" git_log_txt)

  if (NOT (out STREQUAL git_log_txt))
    file(WRITE "${SCRIPT_BIN_DIR}/gitVersion/git_log.txt" out)
    set(LOG_CHANGED TRUE)
  endif()

else()
    file(WRITE "${SCRIPT_BIN_DIR}/gitVersion/git_log.txt" out)
    set(LOG_CHANGED TRUE)
endif()


# Check git status
execute_process(
  COMMAND
  git --git-dir ${SCRIPT_SRC_DIR}/.git --work-tree ${SCRIPT_SRC_DIR}
  status --porcelain applications src testing
  WORKING_DIRECTORY "${SCRIPT_BIN_DIR}"
  RESULT_VARIABLE res
  OUTPUT_VARIABLE out
  )

set(GIT_STATUS "${out}")
string(REPLACE "\n" "\\n" GIT_STATUS "${GIT_STATUS}")

if(GIT_STATUS)
  message(WARNING "Working tree is dirty. Run git status for details.")  
endif()

set(STATUS_CHANGED FALSE)

if (EXISTS "${SCRIPT_BIN_DIR}/gitVersion/git_status.txt")
  file(READ "${SCRIPT_BIN_DIR}/gitVersion/git_status.txt" git_status_txt)

  if (NOT (out STREQUAL git_status_txt))
    file(WRITE "${SCRIPT_BIN_DIR}/gitVersion/git_status.txt" out)
    set(STATUS_CHANGED TRUE)
  endif()

else()
    file(WRITE "${SCRIPT_BIN_DIR}/gitVersion/git_status.txt" out)
    set(STATUS_CHANGED TRUE)

endif()

# Reconfigure gitVersion.cpp if something has changed.
if (LOG_CHANGED OR STATUS_CHANGED)
  configure_file("${SCRIPT_SRC_DIR}/gitVersion/gitVersion.cpp.in" "${SCRIPT_BIN_DIR}/gitVersion/gitVersion.cpp" @ONLY)
endif()
