# test/unit/phys/dca_step/cluster_solver/ctaux/structs

dca_add_gtest(read_write_config_test
        FAST
        GTEST_MPI_MAIN
        INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
        LIBS     FFTW::Double ${DCA_LIBS}
        )
