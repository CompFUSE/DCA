# // Copyright (C) 2022 ETH Zurich
# // Copyright (C) 2022 UT-Battelle, LLC
# // All rights reserved.
# //
# // See LICENSE for terms of usage.
# // See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
# //
# // Author: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
# //
# // This file includes  compiler specific configurations
# // Adapted from QMCPACK/CMake/inspectCompiler.cmake
# // Licensed under the University of Illinois/NCSA Open Source License.
# // Copyright (c) 2022 QMCPACK developers.

# Map the compiler to the internal compiler flag/option customization
# COMPILER is defined upon return
#
# Note: many vendor compilers can be just rebranded customized Clang compiler.
# It requires more recent CMake to handle it properly. We need to handle such cases for older CMake.

# Check compiler version
if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7.0)
  message(STATUS "Compiler Version ${CMAKE_CXX_COMPILER_VERSION}")
  message(FATAL_ERROR "Requires clang 7.0 or higher ")
endif()

if(CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 11.0.0
   AND QMC_CXX_STANDARD EQUAL 17
   AND BUILD_AFQMC)
  message(FATAL_ERROR "Avoid Clang 11.0.0 which cannot compile AFQMC properly with C++17!")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fstrict-aliasing")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fstrict-aliasing")

# treat VLA as error
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror=vla")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wvla")

# set compiler warnings
string(APPEND CMAKE_CXX_FLAGS
       " -Wall -Wno-unused-variable -Wno-overloaded-virtual -Wno-unused-private-field -Wno-unused-local-typedef")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 10.0)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wmisleading-indentation")
endif()

# Set extra optimization specific flags
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffast-math")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")

# Set extra debug flags
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer")

message(STATUS "DCA adds -fstandalone-debug for Debug builds")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fstandalone-debug")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fstandalone-debug")

#--------------------------------------
# Special architectural flags
#--------------------------------------
# case arch
#     x86_64: -march
#     powerpc: -mpcu
#     arm: -mpcu
#     default or cray: none
#--------------------------------------
if(CMAKE_SYSTEM_NAME STREQUAL "CrayLinuxEnvironment")
  # It's a cray machine. Don't do anything
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
  # the case for x86_64
  # check if the user has already specified -march=XXXX option for cross-compiling.
  if(CMAKE_CXX_FLAGS MATCHES "-march=" OR CMAKE_C_FLAGS MATCHES "-march=")
    # make sure that the user specifies -march= for both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS.
    if(CMAKE_CXX_FLAGS MATCHES "-march=" AND CMAKE_C_FLAGS MATCHES "-march=")

    else() #(CMAKE_CXX_FLAGS MATCHES "-march=" AND CMAKE_C_FLAGS MATCHES "-march=")
      message(
        FATAL_ERROR
          "if -march=ARCH is specified by the user, it should be added in both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS!")
    endif() #(CMAKE_CXX_FLAGS MATCHES "-march=" AND CMAKE_C_FLAGS MATCHES "-march=")
  else() #(CMAKE_CXX_FLAGS MATCHES "-march=" OR CMAKE_C_FLAGS MATCHES "-march=")
    # use -march=native
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  endif() #(CMAKE_CXX_FLAGS MATCHES "-march=" OR CMAKE_C_FLAGS MATCHES "-march=")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "ppc64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64")
  # the case for PowerPC and ARM
  # check if the user has already specified -mcpu=XXXX option for cross-compiling.
  if(CMAKE_CXX_FLAGS MATCHES "-mcpu=" OR CMAKE_C_FLAGS MATCHES "-mcpu=")
    # make sure that the user specifies -mcpu= for both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS.
    if(CMAKE_CXX_FLAGS MATCHES "-mcpu=" AND CMAKE_C_FLAGS MATCHES "-mcpu=")

    else() #(CMAKE_CXX_FLAGS MATCHES "-mcpu=" AND CMAKE_C_FLAGS MATCHES "-mcpu=")
      message(
        FATAL_ERROR
          "if -mcpu=ARCH is specified by the user, it should be added in both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS!")
    endif() #(CMAKE_CXX_FLAGS MATCHES "-mcpu=" AND CMAKE_C_FLAGS MATCHES "-mcpu=")
  else() #(CMAKE_CXX_FLAGS MATCHES "-mcpu=" OR CMAKE_C_FLAGS MATCHES "-mcpu=")
    # use -mcpu=native
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mcpu=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=native")
  endif() #(CMAKE_CXX_FLAGS MATCHES "-mcpu=" OR CMAKE_C_FLAGS MATCHES "-mcpu=")
endif()

# Coverage
if(ENABLE_GCOV)
  set(GCOV_SUPPORTED TRUE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
endif(ENABLE_GCOV)

set(XRAY_PROFILE
    FALSE
    CACHE BOOL "Use llvm xray profiling")
set(XRAY_INSTRUCTION_THRESHOLD
    200
    CACHE STRING "Instruction threshold for xray instrumentation")

if(XRAY_PROFILE)
  set(XRAY_FLAGS "-fxray-instrument -fxray-instruction-threshold=${XRAY_INSTRUCTION_THRESHOLD}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${XRAY_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${XRAY_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${XRAY_FLAGS}")
endif(XRAY_PROFILE)
