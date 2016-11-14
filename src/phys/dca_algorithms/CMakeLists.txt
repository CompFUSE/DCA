# DCA(+) algorithms

if (DCA_HAVE_MULDEFS)
  add_library(dca_algorithms STATIC compute_band_structure.cpp)
else()
  add_library(dca_algorithms SHARED compute_band_structure.cpp)
  target_link_libraries(dca_algorithms PRIVATE ${LAPACK_LIBRARIES})
endif()
  
target_include_directories(dca_algorithms PUBLIC ${FFTW_INCLUDE_DIR})
