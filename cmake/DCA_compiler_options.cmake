# Compiler options and tweaks
#
# TODO: - Use target_compile_options().
#       - Use set_property(TARGET target PROPERTY CXX_STANDARD 11) instead of -std=c++11.

# Check for --allow-multiple-definition linker flag.
# References: https://cmake.org/pipermail/cmake/2011-July/045525.html
#             http://projectsymphony.blogspot.ch/2013/03/cmake-linker-test-flags.html
include(CheckCXXCompilerFlag)
set(CMAKE_REQUIRED_FLAGS "-Wl,--allow-multiple-definition")
check_cxx_compiler_flag("" DCA_HAVE_MULDEFS) # Pass empty string since we are testing a linker flag.
if(DCA_HAVE_MULDEFS)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--allow-multiple-definition")
endif()
unset(CMAKE_REQUIRED_FLAGS)

set(WARNINGS "-Wall -Wextra -Wpedantic -Wno-sign-compare")
set(FLAGS "-std=c++14")  # -funroll-loops -finline-functions
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARNINGS} ${FLAGS}")
