
dca_add_gtest(coarsegraining_test
        EXTENSIVE
        GTEST_MAIN
        MPI MPI_NUMPROC 2
        INCLUDE_DIRS ${DCA_INCLUDE_DIRS}
        LIBS         ${DCA_LIBS}
        )
