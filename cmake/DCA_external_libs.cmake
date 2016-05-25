################################################################################
# External libraries
#
# TODO: - Write FindNFFT.cmake.
#       - Write FindSPGLIB.cmake.
#       - Write FindFFTW.cmake.
#       - Use static libraries for NFFT and SPGLIB?
#       - Set DCA_HAVE_XXX to true after XXX was found?
################################################################################

# NFFT
find_library(NFFT_LIBRARY
  NAMES libnfft3.a nfft3
  PATHS ${NFFT_DIR}/lib
  NO_DEFAULT_PATH)

# SPGLIB
find_library(SPGLIB_LIBRARY
  NAMES libsymspg.a symspg
  PATHS ${SPGLIB_DIR}/lib
  NO_DEFAULT_PATH)

# Lapack
if (NOT DCA_LAPACK_IMPLICIT)
  find_package(LAPACK REQUIRED)
endif()

# HDF5
if (NOT DCA_HDF5_IMPLICIT)
  # Find HDF5 by looking for a CMake configuration file (hdf5-1.10.x).
  find_package(HDF5 COMPONENTS C CXX NO_MODULE QUIET)
  if (NOT HDF5_FOUND) 
    # Fall back to a search for a FindHDF5.cmake file and execute it.
    find_package(HDF5 REQUIRED COMPONENTS C CXX)
  endif()
  include_directories(${HDF5_INCLUDE_DIR} ${HDF5_INCLUDE_DIR_CPP})
endif()

# FFTW
if (NOT DCA_FFTW_IMPLICIT)
  find_library(FFTW_LIBRARY NAMES fftw3)
  get_filename_component(FFTW_LIB_DIR ${FFTW_LIBRARY} DIRECTORY)
  get_filename_component(FFTW_DIR     ${FFTW_LIB_DIR} DIRECTORY)
  set(FFTW_INCLUDE_DIR "${FFTW_DIR}/include" CACHE FILEPATH "Path to fftw3.h.")
endif()

set(DCA_EXTERNAL_LIBS
  ${LAPACK_LIBRARIES}
  ${HDF5_EXPORT_LIBRARIES}
  ${HDF5_CXX_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${NFFT_LIBRARY}
  ${FFTW_LIBRARY}
  ${SPRNG_LIBRARY}
)

set(DCA_EXTERNAL_INCLUDES
  ${NFFT_DIR}/include
  ${SPGLIB_DIR}/include
  ${FFTW_INCLUDE_DIR}
  ${HDF5_INCLUDE_DIRS}
  )

mark_as_advanced(
  MPI_LIBRARY MPI_EXTRA_LIBRARY
  NFFT_LIBRARY
  SPGLIB_LIBRARY
  FFTW_INCLUDE_DIR FFTW_LIBRARY
  HDF5_DIR
)

# SPRNG
# Only try to find SPRNG if it is the requested rng to use.
if (${DCA_RNG} STREQUAL "SPRNG")
  # INTERNAL: Is there a find_package for SPRNG?
  dca_add_config_define(DCA_HAVE_SPRNG)

  find_library(SPRNG_LIBRARY
    NAMES libsprng.a sprng
    PATHS ${SPRNG_DIR}/lib
    NO_DEFAULT_PATH)

  if (${SPRNG_LIBRARY} STREQUAL "SPRNG_LIBRARY-NOTFOUND")
    unset(SPRNG_LIBRARY CACHE)
    message(FATAL_ERROR "SPRNG library was not found!\nChoose a different option for the random number generator.")
  endif()
  
  list(APPEND DCA_EXTERNAL_LIBS ${SPRNG_LIBRARY})
  list(APPEND DCA_EXTERNAL_INCLUDES ${SPRNG_DIR}/include)

  mark_as_advanced(SPRNG_LIBRARY)
endif()
