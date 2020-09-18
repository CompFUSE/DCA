################################################################################
# test/integration/cluster-solver/square_lattice/
# CMakeLists.txt
################################################################################
set(TEST_INCLUDES ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR})
set(TEST_LIBS     ${DCA_LIBS};statistical_testing)

dca_add_gtest(ctaux_twoband_Cu_verification_stattest
    STOCHASTIC
    # Run with more ranks for better error detection.
    MPI MPI_NUMPROC 1
    INCLUDE_DIRS ${TEST_INCLUDES}
    LIBS     ${TEST_LIBS}
)

dca_add_gtest(ctaux_twoband_Cu_validation_stattest
        STOCHASTIC
        # Run with more ranks for better error detection.
        MPI MPI_NUMPROC 32
        INCLUDE_DIRS ${TEST_INCLUDES}
        LIBS     ${TEST_LIBS}
)

dca_add_gtest(ctint_twoband_Cu_verification_stattest
        STOCHASTIC
        # Run with more ranks for better error detection.
        MPI MPI_NUMPROC 1
        INCLUDE_DIRS ${TEST_INCLUDES}
        LIBS     ${TEST_LIBS}
)

dca_add_gtest(ctint_twoband_Cu_validation_stattest
        STOCHASTIC
        # Run with more ranks for better error detection.
        MPI MPI_NUMPROC 32
        INCLUDE_DIRS ${TEST_INCLUDES}
        LIBS     ${TEST_LIBS}
)
