
dca_add_gtest(coarsegraining_test
        EXTENSIVE
        MPI MPI_NUMPROC 2
        INCLUDE_DIRS ${DCA_INCLUDE_DIRS}
        LIBS         ${DCA_LIBS}
        )
