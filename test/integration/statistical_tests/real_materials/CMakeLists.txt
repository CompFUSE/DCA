set(TEST_INCLUDES ${DCA_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR})
set(TEST_LIBRARIES ${DCA_LIBS} ${DCA_CUDA_LIBS} statistical_testing)

#  This application was used to build the reference data for the NiO test.
#add_executable(NiO_fullDCA NiO_fullDCA.cpp)
#target_link_libraries(NiO_fullDCA ${TEST_LIBRARIES})
#target_include_directories(NiO_fullDCA PRIVATE ${TEST_INCLUDES})
#target_compile_definitions(NiO_fullDCA PRIVATE TEST_DIRECTORY="${CMAKE_CURRENT_SOURCE_DIR}/")

dca_add_gtest(NiO_stattest
    EXTENSIVE
    PTHREADS
    # Run with more ranks for better error detection.
    MPI MPI_NUMPROC 1
    INCLUDE_DIRS "${TEST_INCLUDES}"
    LIBS     "${TEST_LIBRARIES}"
    )

configure_file(NiO_coarse_grained.hdf5 ${CMAKE_CURRENT_BINARY_DIR}/NiO_coarse_grained.hdf5 COPYONLY)
configure_file(NiO_covariance_input.hdf5 ${CMAKE_CURRENT_BINARY_DIR}/NiO_covariance_input.hdf5 COPYONLY)
