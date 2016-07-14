# gaussian quadrature

# We need to create a shared library to avoid a multiple definition error when linking.
add_library(gaussian_quadrature SHARED compute_weights_and_abscissas.cpp simplex_gm_rule.cpp)

# simplex_gm_rule.cpp includes simplex_gm_rule.hpp without any path prefix. We therefore need to add
# the full path to the include directories.
target_include_directories(gaussian_quadrature
                           PRIVATE ${PROJECT_SOURCE_DIR}/include
                           PRIVATE ${PROJECT_SOURCE_DIR}/include/dca/math/geometry/gaussian_quadrature
                           PRIVATE ${FFTW_INCLUDE_DIR})

# The shared library already needs to be linked to LAPACK.
target_link_libraries(gaussian_quadrature PRIVATE ${LAPACK_LIBRARIES})

# Silence compiler warnings caused by simplex_gm_rule.cpp.
target_compile_options(gaussian_quadrature PRIVATE "-Wno-unused-parameter")
if (CMAKE_COMPILER_IS_GNUCXX)
target_compile_options(gaussian_quadrature PRIVATE "-Wno-unused-but-set-variable")
endif()
