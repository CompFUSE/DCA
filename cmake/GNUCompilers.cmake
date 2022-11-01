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

# Check compiler version
if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0)
  message(FATAL_ERROR "Requires gcc 8.0 or higher ")
endif()

# Set gnu specific flags (which we always want)
add_definitions(-Drestrict=__restrict__)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -finline-limit=1000 -fstrict-aliasing -funroll-all-loops")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -finline-limit=1000 -fstrict-aliasing -funroll-all-loops")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer")

# Suppress compile warnings
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-deprecated")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")

# treat VLA as error
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror=vla")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wvla")

# set compiler warnings
set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -Wcomment -Wmisleading-indentation -Wmaybe-uninitialized -Wuninitialized -Wreorder -Wno-unknown-pragmas -Wno-sign-compare"
)

# Set extra optimization specific flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffast-math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")

