# File: cmake/gitVersion_rebuild.cmake
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# This cmake script is executed every time 'make' is executed.
# It checks the last git log (git log -1) and git status (git status --porcelain).
# If either of them changes, gitVersion.cpp will be reconfigured which triggers the
# recompilation of libgitVersion.

# Check git log
execute_process(
  COMMAND
  git
  log -1
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"  # .../build
  RESULT_VARIABLE res
  OUTPUT_VARIABLE out
  )

set(GIT_LOG "${out}")
# Replace newlines with characters "\n" for use in std::string.
string(REPLACE "\n" "\\n" GIT_LOG "${GIT_LOG}")

set(LOG_CHANGED FALSE)

if (EXISTS "${CMAKE_BINARY_DIR}/gitVersion/git_log.txt")
  file(READ "${CMAKE_BINARY_DIR}/gitVersion/git_log.txt" git_log_txt)

  if (NOT (out STREQUAL git_log_txt))
    file(WRITE "${CMAKE_BINARY_DIR}/gitVersion/git_log.txt" out)
    set(LOG_CHANGED TRUE)
  endif()

else()
    file(WRITE "${CMAKE_BINARY_DIR}/gitVersion/git_log.txt" out)
    set(LOG_CHANGED TRUE)
endif()


# Check git status
execute_process(
  COMMAND
  git
  status --porcelain
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
  RESULT_VARIABLE res
  OUTPUT_VARIABLE out
  )

set(GIT_STATUS "${out}")
string(REPLACE "\n" "\\n" GIT_STATUS "${GIT_STATUS}")

if(GIT_STATUS)
  message(WARNING "Working tree is dirty. Run git status for details.")  
endif()

set(STATUS_CHANGED FALSE)

if (EXISTS "${CMAKE_BINARY_DIR}/gitVersion/git_status.txt")
  file(READ "${CMAKE_BINARY_DIR}/gitVersion/git_status.txt" git_status_txt)

  if (NOT (out STREQUAL git_status_txt))
    file(WRITE "${CMAKE_BINARY_DIR}/gitVersion/git_status.txt" out)
    set(STATUS_CHANGED TRUE)
  endif()

else()
    file(WRITE "${CMAKE_BINARY_DIR}/gitVersion/git_status.txt" out)
    set(STATUS_CHANGED TRUE)

endif()

# Reconfigure gitVersion.cpp if something has changed.
if (LOG_CHANGED OR STATUS_CHANGED)
  configure_file("${CMAKE_SOURCE_DIR}/gitVersion/gitVersion.cpp.in" "${CMAKE_BINARY_DIR}/gitVersion/gitVersion.cpp" @ONLY)
endif()
