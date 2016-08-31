################################################################################
# Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
#
# Checks for external libraries and creates a global lists of them and the corresponding include paths.
#
# TODO: - Set DCA_HAVE_XXX to true after XXX was found?
#       - Use static or shared libraries?

set(DCA_EXTERNAL_LIBS "" CACHE INTERNAL "")
set(DCA_EXTERNAL_INCLUDES "" CACHE INTERNAL "")

################################################################################
# NFFT
find_library(NFFT_LIBRARY nfft3 HINTS ${NFFT_DIR}/lib)
find_path(NFFT_INCLUDE_DIR nfft3.h HINTS ${NFFT_DIR}/include)

mark_as_advanced(NFFT_LIBRARY NFFT_INCLUDE_DIR)

list(APPEND DCA_EXTERNAL_LIBS ${NFFT_LIBRARY})
list(APPEND DCA_EXTERNAL_INCLUDES ${NFFT_INCLUDE_DIR})

################################################################################
# SPGLIB
find_library(SPGLIB_LIBRARY symspg HINTS ${SPGLIB_DIR}/lib)
find_path(SPGLIB_INCLUDE_DIR spglib.h HINTS ${SPGLIB_DIR}/include/spglib)

mark_as_advanced(SPGLIB_LIBRARY SPGLIB_INCLUDE_DIR)

list(APPEND DCA_EXTERNAL_LIBS ${SPGLIB_LIBRARY})
list(APPEND DCA_EXTERNAL_INCLUDES ${SPGLIB_INCLUDE_DIR})

################################################################################
# Lapack
if (NOT DCA_HAVE_LAPACK)
  find_package(LAPACK REQUIRED)
endif()

mark_as_advanced(LAPACK_LIBRARIES)

list(APPEND DCA_EXTERNAL_LIBS ${LAPACK_LIBRARIES})

################################################################################
# HDF5
if (NOT DCA_HAVE_HDF5)
  # Find HDF5 by looking for a CMake configuration file (hdf5-1.10.x).
  find_package(HDF5 COMPONENTS C CXX NO_MODULE QUIET)

  if (NOT HDF5_FOUND)
    # Fall back to a search for a FindHDF5.cmake file and execute it.
    find_package(HDF5 REQUIRED COMPONENTS C CXX)
  endif()
endif()

mark_as_advanced(HDF5_CXX_INCLUDE_DIR HDF5_C_INCLUDE_DIR HDF5_DIR)

list(APPEND DCA_EXTERNAL_LIBS ${HDF5_LIBRARIES})
list(APPEND DCA_EXTERNAL_INCLUDES ${HDF5_INCLUDE_DIRS})

################################################################################
# FFTW
set(FFTW_INCLUDE_DIR "" CACHE PATH "Path to fftw3.h.")
set(FFTW_LIBRARY "" CACHE FILEPATH "Path to FFTW3 library.")
mark_as_advanced(FFTW_INCLUDE_DIR FFTW_LIBRARY)

if (NOT DCA_HAVE_FFTW)
  if (NOT FFTW_INCLUDE_DIR OR NOT FFTW_LIBRARY)
    message(FATAL_ERROR "FFTW_INCLUDE_DIR and FFTW_LIBRARY have to be set.")
  endif()
endif()

list(APPEND DCA_EXTERNAL_LIBS ${FFTW_LIBRARY})
list(APPEND DCA_EXTERNAL_INCLUDES ${FFTW_INCLUDE_DIR})

################################################################################
# SPRNG
if (DCA_RNG MATCHES "^SPRNG")
  if (NOT DCA_HAVE_SPRNG)
    message(FATAL_ERROR "SPRNG library was not found! Choose a different random number generator.")
  endif()

  list(APPEND DCA_EXTERNAL_LIBS ${SPRNG_LIBRARY})
  list(APPEND DCA_EXTERNAL_INCLUDES ${SPRNG_INCLUDE_DIR})
endif()

################################################################################
# message("DCA_EXTERNAL_LIBS = ${DCA_EXTERNAL_LIBS}")
# message("DCA_EXTERNAL_INCLUDES = ${DCA_EXTERNAL_INCLUDES}")
