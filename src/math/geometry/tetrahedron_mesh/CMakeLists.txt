# Tetrahedron mesh

if (DCA_HAVE_MULDEFS)
  add_library(tetrahedron_mesh STATIC facet.cpp tetrahedron.cpp)
else()
  add_library(tetrahedron_mesh SHARED facet.cpp tetrahedron.cpp)
  target_link_libraries(tetrahedron_mesh PRIVATE ${GNUPLOT_INTERFACE_LIBRARY} ${LAPACK_LIBRARIES})
endif()
  target_include_directories(tetrahedron_mesh PUBLIC ${FFTW_INCLUDE_DIR})
