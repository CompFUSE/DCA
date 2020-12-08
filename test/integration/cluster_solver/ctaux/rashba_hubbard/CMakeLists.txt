dca_add_gtest(ctaux_rashba_hubbard_test
        EXTENSIVE
        CUDA
        INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
        LIBS     ${DCA_LIBS}
        )
