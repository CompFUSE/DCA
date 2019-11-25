################################################################################
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# Compiler options and tweaks
#
# TODO: - Use target_compile_options().
#       - Use set_property(TARGET target PROPERTY CXX_STANDARD 14) instead of -std=c++14.

# Check for --allow-multiple-definition linker flag.
# References: https://cmake.org/pipermail/cmake/2011-July/045525.html
#             http://projectsymphony.blogspot.ch/2013/03/cmake-linker-test-flags.html
# include(CheckCXXCompilerFlag)
# set(CMAKE_REQUIRED_FLAGS "-Wl,--allow-multiple-definition")
# check_cxx_compiler_flag("" DCA_HAVE_MULDEFS) # Pass empty string since we are testing a linker flag.
# if (DCA_HAVE_MULDEFS)
#   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--allow-multiple-definition")
# endif()
# unset(CMAKE_REQUIRED_FLAGS)

# Warnings
set(DCA_WARNINGS -Wall -Wextra -Wpedantic -Wno-sign-compare -Wno-dangling-else)

# Languange standard
set(DCA_STD_FLAG -std=c++17)

# Set C and CXX flags.
add_compile_options(${DCA_WARNINGS})
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:${DCA_STD_FLAG}>")

# Add support for multithreading with the Pthread library.
add_compile_options("-pthread")

# Set NVCC flags.
if (DCA_HAVE_CUDA)
  if (CUDA_VERSION VERSION_GREATER "8.1.0")
    list(APPEND CUDA_NVCC_FLAGS
      -arch=${CUDA_GPU_ARCH}
      -std=c++14
      -Xcompiler -Wall
      -Xcompiler -Wextra
      -Xcompiler -Wno-unused-parameter
      -Xcompiler -Wno-switch
      -Xcompiler ${DCA_THREADING_FLAGS})
  else (CUDA_VERSION VERSION_GREATER "8.1.0")
    list(APPEND CUDA_NVCC_FLAGS
      -arch=${CUDA_GPU_ARCH}
      -std=c++11
      -Xcompiler -Wall
      -Xcompiler -Wextra
      -Xcompiler -Wno-unused-parameter
      -Xcompiler -Wno-switch
      -Xcompiler ${DCA_THREADING_FLAGS})
  endif (CUDA_VERSION VERSION_GREATER "8.1.0")
endif()
