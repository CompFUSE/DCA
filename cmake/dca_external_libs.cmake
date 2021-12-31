################################################################################
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# Checks for external libraries and creates a global lists of them and the corresponding include
# paths.
# We explicitely search for static libraries because compute nodes cannot load shared libraries from
# e.g. the project directory.

set(DCA_EXTERNAL_LIBS "" CACHE INTERNAL "")
set(DCA_EXTERNAL_INCLUDE_DIRS "" CACHE INTERNAL "")

################################################################################
# Lapack
if (NOT DCA_HAVE_LAPACK)
  find_package(LAPACK REQUIRED)
else()
  add_library(LAPACK::LAPACK INTERFACE IMPORTED)
  set_target_properties(LAPACK::LAPACK PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${LAPACK_INCLUDE_DIRS}"
  INTERFACE_COMPILE_DEFINITION "${LAPACK_DEFINITIONS}"
  IMPORTED_LINK_INTERFACE_LANGUAGES "C;CXX"
  IMPORTED_LOCATION "${LAPACK_LIBRARY}")
endif()

mark_as_advanced(LAPACK_LIBRARIES)
message("LAPACK_INCLUDE_DIRS: ${LAPACK_INCLUDE_DIRS}")
message("LAPACK_LIBRARIES: ${LAPACK_FOUND} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES}")
list(APPEND DCA_EXTERNAL_LIBS ${LAPACK_LIBRARIES})

################################################################################
# Blas
if (NOT DCA_HAVE_BLAS)
  find_package(BLAS REQUIRED)
endif()

mark_as_advanced(LAPACK_LIBRARIES)
message("LAPACK_LIBRARIES: ${LAPACK_FOUND} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES}")
list(APPEND DCA_EXTERNAL_LIBS ${LAPACK_LIBRARIES})

################################################################################
# HDF5

find_package(HDF5 REQUIRED COMPONENTS C CXX)

################################################################################
# ADIOS2
if (DCA_WITH_ADIOS2)
set(DCA_HAVE_ADIOS2 FALSE CACHE INTERNAL "")
find_package(ADIOS2)
if (ADIOS2_FOUND)
  set(DCA_HAVE_ADIOS2 TRUE CACHE INTERNAL "")
  #message("ADIOS2: libraries ${ADIOS2_LIBRARIES}")
endif()
endif()

################################################################################
################################################################################
# FFTW
find_package(FFTW REQUIRED)

list(APPEND DCA_EXTERNAL_LIBS ${FFTW_LIBRARY})
list(APPEND DCA_EXTERNAL_INCLUDE_DIRS ${FFTW_INCLUDE_DIR})

################################################################################
# Simplex GM Rule
add_subdirectory(${PROJECT_SOURCE_DIR}/libs/simplex_gm_rule)
# list(APPEND DCA_EXTERNAL_LIBS ${SIMPLEX_GM_RULE_LIBRARY})
# list(APPEND DCA_EXTERNAL_INCLUDE_DIRS ${SIMPLEX_GM_RULE_INCLUDE_DIR})

################################################################################
# Gnuplot
list(APPEND DCA_EXTERNAL_LIBS ${GNUPLOT_INTERFACE_LIBRARY})
list(APPEND DCA_EXTERNAL_INCLUDE_DIRS ${GNUPLOT_INTERFACE_INCLUDE_DIR})

message("DCA_EXTERNAL_LIBS = ${DCA_EXTERNAL_LIBS}")
message("DCA_EXTERNAL_INCLUDE_DIRS = ${DCA_EXTERNAL_INCLUDE_DIRS}")
