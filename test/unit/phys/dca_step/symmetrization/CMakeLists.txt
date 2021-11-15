# symmetrization

dca_add_gtest(symmetrize_test
        FAST
        GTEST_MAIN
        INCLUDE_DIRS ${DCA_INCLUDES};${PROJECT_SOURCE_DIR};${FFTW_INCLUDE_DIR}
        LIBS  ${DCA_LIBS}
        )
