# Coarse-graining

if (DCA_HAVE_MULDEFS)
  add_library(coarsegraining
    STATIC
    coarsegrain_domain_names.cpp
    tetrahedron_routines_harmonic_function.cpp)

else()
  add_library(coarsegraining
    SHARED
    coarsegrain_domain_names.cpp
    tetrahedron_routines_harmonic_function.cpp)
  target_link_libraries(coarsegraining PUBLIC tetrahedron_mesh)
endif()
  
target_include_directories(coarsegraining PUBLIC ${FFTW_INCLUDE_DIR})
