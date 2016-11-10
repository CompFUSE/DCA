# Tetrahedron mesh

add_library(tetrahedron_mesh STATIC facet.cpp tetrahedron.cpp)
target_include_directories(tetrahedron_mesh PUBLIC ${FFTW_INCLUDE_DIR})
