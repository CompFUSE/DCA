if (DCA_THREADING_LIBRARY STREQUAL POSIX OR
    DCA_THREADING_LIBRARY STREQUAL stdthread)
  add_subdirectory(thread_qmci)
endif()

add_subdirectory(util)
