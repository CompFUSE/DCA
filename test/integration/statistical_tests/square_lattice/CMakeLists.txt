################################################################################
# test/integration/cluster-solver/square_lattice/
# CMakeLists.txt
################################################################################
set(TEST_INCLUDES ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR})
set(TEST_LIBS     ${DCA_LIBS};statistical_testing)

dca_add_gtest(ctaux_square_lattice_verification_stattest
    STOCHASTIC
    # Run with more ranks for better error detection.
    MPI MPI_NUMPROC 1
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${TEST_LIBS}
    )

dca_add_gtest(ctaux_square_lattice_validation_stattest
    STOCHASTIC
    # Run with more ranks for better error detection.
    MPI MPI_NUMPROC 32
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${TEST_LIBS}
    )

dca_add_gtest(ctint_square_lattice_verification_stattest
    STOCHASTIC
    # Run with more ranks for better error detection.
    MPI MPI_NUMPROC 1
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${TEST_LIBS}
    )

dca_add_gtest(ctint_square_lattice_validation_stattest
    STOCHASTIC
    # Run with more ranks for better error detection.
    MPI MPI_NUMPROC 32
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${TEST_LIBS}
    )


configure_file(data.ed.hdf5 ${CMAKE_CURRENT_BINARY_DIR}/data.ed.hdf5 COPYONLY)
configure_file(verification_covariance_input.hdf5 ${CMAKE_CURRENT_BINARY_DIR}/verification_covariance_input.hdf5 COPYONLY)
