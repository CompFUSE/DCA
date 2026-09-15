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

# Disorder clean-limit stat test (Step 1): drives the r-space disorder pipeline at V(R) = 0 and
# checks it reproduces the clean ED reference. Only built with the disorder option, never for GPU.
if(DCA_WITH_DISORDER AND NOT DCA_WITH_CUDA AND NOT DCA_WITH_HIP)
  dca_add_gtest(disorder_ctaux_square_lattice_stattest
      STOCHASTIC
      # Use the machine's core count rather than a fixed number; more ranks give better statistics.
      MPI MPI_NUMPROC ${CPUS}
      INCLUDE_DIRS ${TEST_INCLUDES}
      LIBS     ${TEST_LIBS}
      )

  # Disorder finite-V robustness test (Step 2): runs the binary + box disorder ensemble at V=0.1 and
  # checks G_k_w stays finite/sane. A smoke test (no accuracy/baseline); grouped with the stochastic
  # tests so it is not built unless DCA_WITH_TESTS_STOCHASTIC is on.
  dca_add_gtest(disorder_finite_v_stattest
      STOCHASTIC
      MPI MPI_NUMPROC 1
      INCLUDE_DIRS ${TEST_INCLUDES}
      LIBS     ${TEST_LIBS}
      )
endif()

configure_file(data.ed.hdf5 ${CMAKE_CURRENT_BINARY_DIR}/data.ed.hdf5 COPYONLY)
configure_file(verification_covariance_input.hdf5 ${CMAKE_CURRENT_BINARY_DIR}/verification_covariance_input.hdf5 COPYONLY)
