################################################################################
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# Checks for external libraries and creates a global lists of them and the corresponding include
# paths.
# We explicitely search for static libraries because compute nodes cannot load shared libraries from
# e.g. the project directory.
#
#
# TODO: - Set DCA_HAVE_XXX to true after XXX was found?
#       - Use static or shared libraries?

set(DCA_EXTERNAL_LIBS "" CACHE INTERNAL "")
set(DCA_EXTERNAL_INCLUDE_DIRS "" CACHE INTERNAL "")

################################################################################
# Lapack
if (NOT DCA_HAVE_LAPACK)
  find_package(LAPACK REQUIRED)
endif()

mark_as_advanced(LAPACK_LIBRARIES)

list(APPEND DCA_EXTERNAL_LIBS ${LAPACK_LIBRARIES})

################################################################################
# HDF5
set(HDF5_ROOT "" CACHE PATH "Path to HDF5 installation directory.")

if (NOT DCA_HAVE_HDF5)
  # Find HDF5 by looking for a CMake configuration file (hdf5-1.10.x).
  find_package(HDF5 COMPONENTS C CXX NO_MODULE QUIET)

  if (NOT HDF5_FOUND)
    # Fall back to a search for a FindHDF5.cmake file and execute it.
    find_package(HDF5 REQUIRED COMPONENTS C CXX)
  endif()
endif()

mark_as_advanced(HDF5_ROOT HDF5_CXX_INCLUDE_DIR HDF5_C_INCLUDE_DIR HDF5_DIR)

list(APPEND DCA_EXTERNAL_LIBS ${HDF5_LIBRARIES})
list(APPEND DCA_EXTERNAL_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})

################################################################################
# FFTW
set(FFTW_INCLUDE_DIR "" CACHE PATH "Path to fftw3.h.")
set(FFTW_LIBRARY "" CACHE FILEPATH "FFTW3 library.")
mark_as_advanced(FFTW_INCLUDE_DIR FFTW_LIBRARY)

if (NOT DCA_HAVE_FFTW)
  if (NOT FFTW_INCLUDE_DIR OR NOT FFTW_LIBRARY)
    message(FATAL_ERROR "FFTW_INCLUDE_DIR and FFTW_LIBRARY have to be set.")
  endif()
endif()

list(APPEND DCA_EXTERNAL_LIBS ${FFTW_LIBRARY})
list(APPEND DCA_EXTERNAL_INCLUDE_DIRS ${FFTW_INCLUDE_DIR})

################################################################################
# SPRNG
if (DCA_RNG MATCHES "^SPRNG")
  if (NOT DCA_HAVE_SPRNG)
    message(FATAL_ERROR "SPRNG library was not found! Choose a different random number generator.")
  endif()

  list(APPEND DCA_EXTERNAL_LIBS ${SPRNG_LIBRARY})
  list(APPEND DCA_EXTERNAL_INCLUDE_DIRS ${SPRNG_INCLUDE_DIR})
endif()

################################################################################
# message("DCA_EXTERNAL_LIBS = ${DCA_EXTERNAL_LIBS}")
# message("DCA_EXTERNAL_INCLUDE_DIRS = ${DCA_EXTERNAL_INCLUDE_DIRS}")
