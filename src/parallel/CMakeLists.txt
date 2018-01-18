if (DCA_THREADING_LIBRARY STREQUAL POSIX)
  add_subdirectory(pthreading)
endif()

add_subdirectory(util)
