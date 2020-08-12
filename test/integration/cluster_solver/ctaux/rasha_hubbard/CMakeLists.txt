dca_add_gtest(ctaux_rasha_hubbard_test
        EXTENSIVE
        # MPI MPI_NUMPROC 4
        INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
        LIBS     ${DCA_LIBS}
        )
