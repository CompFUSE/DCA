################################################################################
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# This cmake script is executed every time 'make' is executed. It checks the last git log
# (git log -1) and git status (git status --porcelain applications include src test). If one of the
# two changes, git_version.cpp will be reconfigured, which triggers the recompilation of
# libGitVersion.
#
# SCRIPT_SRC_DIR is passed in from the custom rule as the PROJECT_SOURCE_DIR is not valid when run
# as a rule.

# Check git log
execute_process(
  COMMAND git --git-dir ${SCRIPT_SRC_DIR}/.git --work-tree ${SCRIPT_SRC_DIR} log -1
  WORKING_DIRECTORY "${SCRIPT_SRC_DIR}"
  RESULT_VARIABLE res
  OUTPUT_VARIABLE out)

set(GIT_LOG "${out}")
# Replace newlines with characters "\n" for use in std::string.
string(REPLACE "\n" "\\n" GIT_LOG "${GIT_LOG}")
# Escape double quotes within the string literals.
string(REPLACE "\"" "\\\"" GIT_LOG "${GIT_LOG}")

set(LOG_CHANGED FALSE)

if (EXISTS "${SCRIPT_BIN_DIR}/src/util/git_log.txt")
  file(READ "${SCRIPT_BIN_DIR}/src/util/git_log.txt" git_log_txt)

  if (NOT (out STREQUAL git_log_txt))
    file(WRITE "${SCRIPT_BIN_DIR}/src/util/git_log.txt" "${out}")
    set(LOG_CHANGED TRUE)
  endif()

else()
    file(WRITE "${SCRIPT_BIN_DIR}/src/util/git_log.txt" "${out}")
    set(LOG_CHANGED TRUE)
endif()


# Check git status
execute_process(
  COMMAND git --git-dir ${SCRIPT_SRC_DIR}/.git --work-tree ${SCRIPT_SRC_DIR}
  status --porcelain applications include src test
  WORKING_DIRECTORY "${SCRIPT_SRC_DIR}"
  RESULT_VARIABLE res
  OUTPUT_VARIABLE out)

set(GIT_STATUS "${out}")
string(REPLACE "\n" "\\n" GIT_STATUS "${GIT_STATUS}")
string(REPLACE "\"" "\\\"" GIT_STATUS "${GIT_STATUS}")

if(GIT_STATUS)
  message(WARNING "Working tree is dirty. Run git status for details.")
endif()

set(STATUS_CHANGED FALSE)

if (EXISTS "${SCRIPT_BIN_DIR}/src/util/git_status.txt")
  file(READ "${SCRIPT_BIN_DIR}/src/util/git_status.txt" git_status_txt)

  if (NOT (out STREQUAL git_status_txt))
    file(WRITE "${SCRIPT_BIN_DIR}/src/util/git_status.txt" "${out}")
    set(STATUS_CHANGED TRUE)
  endif()

else()
    file(WRITE "${SCRIPT_BIN_DIR}/src/util/git_status.txt" "${out}")
    set(STATUS_CHANGED TRUE)

endif()

# Reconfigure git_version.cpp if something has changed.
if (LOG_CHANGED OR STATUS_CHANGED)
  configure_file("${SCRIPT_SRC_DIR}/src/util/git_version.cpp.in"
    "${SCRIPT_BIN_DIR}/src/util/git_version.cpp" @ONLY)
endif()
