################################################################################
# External libraries
#
# TODO: - Write FindNFFT.cmake.
#       - Write FindSPGLIB.cmake.
#       - Write FindFFTW.cmake.
#       - Use static libraries for NFFT and SPGLIB?
#       - Set DCA_XXX_AVAILABLE to true after XXX was found?
################################################################################

# NFFT
find_library(NFFT_LIBRARY
  NAMES libnfft3.a nfft3
  PATHS ${EBROOTNFFT}/lib
  NO_DEFAULT_PATH
  )

# SPGLIB
find_library(SPGLIB_LIBRARY
  NAMES libsymspg.a symspg
  PATHS ${EBROOTSPGLIB}/lib
  NO_DEFAULT_PATH
  )

# Lapack
if (NOT DCA_LAPACK_IMPLICIT)
  find_package(LAPACK REQUIRED)
endif()

# # Scalapack
# find_library(VECLIBFORT NAMES veclibFort)
# find_library(SCALAPACK  NAMES scalapack)

# HDF5
if (NOT DCA_HDF5_IMPLICIT)
  find_package(HDF5 REQUIRED COMPONENTS CXX)
  # if(NOT HDF5_FOUND)
  #    set(HDF5_DIR "${CMAKE_SOURCE_DIR}/libs/hdf5-1.8.11")
  #    find_library(HDF5_LIBRARIES      NAMES hdf5      PATHS ${HDF5_DIR}/build/lib NO_DEFAULT_PATH)
  #    find_library(HDF5_CXX_LIBRARIES  NAMES hdf5_cpp  PATHS ${HDF5_DIR}/build/lib NO_DEFAULT_PATH)
  # endif()
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
  ${HDF5_LIBRARIES}
  ${HDF5_CXX_LIBRARIES}
  ${NFFT_LIBRARY}
  ${FFTW_LIBRARY}
  ${SPGLIB_LIBRARY}
  )
# message("DCA_EXTERNAL_LIBS: ${DCA_EXTERNAL_LIBS}")

set(DCA_EXTERNAL_INCLUDES
  ${EBROOTNFFT}/include
  ${EBROOTSPGLIB}/include
  ${FFTW_INCLUDE_DIR}
  ${HDF5_INCLUDE_DIRS}
  )
# message("DCA_EXTERNAL_INCLUDES: ${DCA_EXTERNAL_INCLUDES}")
