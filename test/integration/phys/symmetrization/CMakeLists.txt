# Symmetrization characterization harness

# Test binaries for derived symmetrization. One binary per (lattice, point group)
# case: the DCA domain + cluster-symmetry singletons are keyed on dimension/role,
# not on the lattice, so cases of the same dimension cannot coexist in one process.
# Each case is selected from the shared source via TEST_DEFINES.
dca_add_gtest(symmetrize_characterization_square_D4_test
        FAST
        GTEST_MPI_MAIN
        INCLUDE_DIRS ${DCA_INCLUDES};${PROJECT_SOURCE_DIR};${FFTW_INCLUDE_DIR}
        LIBS  ${DCA_LIBS}
        CUSTOM_SOURCE symmetrize_characterization_test.cpp
        TEST_DEFINES SquareD4
        )

dca_add_gtest(symmetrize_characterization_threeband_D4_test
        FAST
        GTEST_MPI_MAIN
        INCLUDE_DIRS ${DCA_INCLUDES};${PROJECT_SOURCE_DIR};${FFTW_INCLUDE_DIR}
        LIBS  ${DCA_LIBS}
        CUSTOM_SOURCE symmetrize_characterization_test.cpp
        TEST_DEFINES ThreebandD4
        )

# Regression case for deriveOrbitalOpForOp's -1 guard: Kagome produces the -1 which the solver
# must skip instead of indexing H0 out of bounds.
dca_add_gtest(symmetrize_characterization_kagome_test
        FAST
        GTEST_MPI_MAIN
        INCLUDE_DIRS ${DCA_INCLUDES};${PROJECT_SOURCE_DIR};${FFTW_INCLUDE_DIR}
        LIBS  ${DCA_LIBS}
        CUSTOM_SOURCE symmetrize_characterization_test.cpp
        TEST_DEFINES KagomeNoSym
        )

# The motivating case for deriving P from H0: FeAs's two orbitals share a site, so the geometric
# band image is always the identity and half of D4 is recoverable only by searching for the swap.
dca_add_gtest(symmetrize_characterization_fe_as_test
        FAST
        GTEST_MPI_MAIN
        INCLUDE_DIRS ${DCA_INCLUDES};${PROJECT_SOURCE_DIR};${FFTW_INCLUDE_DIR}
        LIBS  ${DCA_LIBS}
        CUSTOM_SOURCE symmetrize_characterization_test.cpp
        TEST_DEFINES FeAs
        )

dca_add_gtest(symmetrize_characterization_singleband_chain_test
        FAST
        GTEST_MPI_MAIN
        INCLUDE_DIRS ${DCA_INCLUDES};${PROJECT_SOURCE_DIR};${FFTW_INCLUDE_DIR}
        LIBS  ${DCA_LIBS}
        CUSTOM_SOURCE symmetrize_characterization_test.cpp
        TEST_DEFINES SinglebandChain
        )
