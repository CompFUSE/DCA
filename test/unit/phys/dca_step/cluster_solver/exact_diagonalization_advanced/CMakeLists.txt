set(TEST_LIBS function;${LAPACK_LIBRARIES};quantum_domains;cluster_domains;
              time_and_frequency_domains;json;dca_algorithms)

dca_add_gtest(hamiltonian_test
    FAST
    GTEST_MAIN
    INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
    LIBS function json time_and_frequency_domains cluster_domains quantum_domains ${LAPACK_LIBRARIES}
    )

dca_add_gtest(non_density_int_hamiltonian_test
    FAST
    GTEST_MAIN
    INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
    LIBS ${TEST_LIBS}
    )
