# test/unit/phys/dca_step/cluster_solver/shared_tools/tp

dca_add_gtest(hdf5_g4_file_io_test
  FAST
  GTEST_MPI_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     ${DCA_LIBS}
  )

if (DCA_HAVE_ADIOS2)
  dca_add_gtest(adios2_g4_file_io_test
    FAST
    INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
    LIBS     ${DCA_LIBS}
  )

endif()

# Baseline write is broken
# dca_add_gtest(tp_accumulator_test
#   FAST
#   GTEST_MPI_MAIN
#   INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
#   LIBS enumerations function timer parallel_stdthread parallel_util dca_domains dca_math modern_string_utils dca_io parallel_no_concurrency
#   )

dca_add_gtest(tp_accumulator_test_complex
  FAST
  GTEST_MPI_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     ${DCA_LIBS}
  )

dca_add_gtest(tp_accumulator_singleband_test
  FAST
  GTEST_MPI_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     ${DCA_LIBS}
  )

set(TPACC_PARTICLE_HOLE_LIBS gpu_utils cluster_domains time_and_frequency_domains quantum_domains timer json random
  enumerations dca_algorithms parallel_stdthread parallel_util ${LAPACK_LIBRARIES})

if (DCA_HAVE_ADIOS2)
  set(TPACC_PARTICLE_HOLE_LIBS dca_adios2 ${TPACC_PARTICLE_HOLE_LIBS})
endif()

# dca_add_gtest(tp_accumulator_particle_hole_test
#   FAST
#   CUDA
#   GTEST_MPI_MAIN
#   INCLUDE_DIRS ${PROJECT_SOURCE_DIR} ${FFTW_INCLUDE_DIR} ${HDF5_INCLUDE_DIR}
#   LIBS function FFTW::Double cluster_domains time_and_frequency_domains quantum_domains timer json random
#        enumerations dca_algorithms ${DCA_THREADING_LIBS} dca_io parallel_util ${LAPACK_LIBRARIES} gpu_utils parallel_no_concurrency)

dca_add_gtest(tp_accumulator_gpu_test
  FAST
  CUDA
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     ${DCA_LIBS} ${DCA_KERNEL_LIBS}
  CUSTOM_SOURCE tp_accumulator_gpu_test.cpp
  TEST_DEFINES DoubleSquare

)

dca_add_gtest(tp_accumulator_rashba_gpu_test
  FAST
  CUDA
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     ${DCA_LIBS} ${DCA_KERNEL_LIBS}
  CUSTOM_SOURCE tp_accumulator_gpu_test.cpp
  TEST_DEFINES  ComplexDoubleRashba
  )

dca_add_gtest(tp_accumulator_moire_hubbard_gpu_test
  FAST
  CUDA
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     ${DCA_LIBS} ${DCA_KERNEL_LIBS}
  CUSTOM_SOURCE tp_accumulator_gpu_test.cpp
  TEST_DEFINES  ComplexDoubleMoireHubbard
  )

dca_add_gtest(tp_accumulator_hund_gpu_test
  FAST
  CUDA
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     ${DCA_LIBS} ${DCA_KERNEL_LIBS}
  CUSTOM_SOURCE tp_accumulator_gpu_test.cpp
  TEST_DEFINES  DoubleHund
  )

dca_add_gtest(tp_accumulator_DoubleKagome_gpu_test
  FAST
  CUDA
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     ${DCA_LIBS} ${DCA_KERNEL_LIBS}
  CUSTOM_SOURCE tp_accumulator_gpu_test.cpp
  TEST_DEFINES  DoubleKagome
  )

dca_add_gtest(tp_accumulator_FeAs_gpu_test
  FAST
  CUDA
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     ${DCA_LIBS} ${DCA_KERNEL_LIBS}
  CUSTOM_SOURCE tp_accumulator_gpu_test.cpp
  TEST_DEFINES  DoubleFeAs
  )

dca_add_gtest(tp_accumulator_complex_g0_gpu_test
  FAST
  CUDA
        INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
        LIBS     ${DCA_LIBS} ${DCA_KERNEL_LIBS}
        )

dca_add_gtest(tp_accumulator_singleband_gpu_test
  FAST
  CUDA
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     ${DCA_LIBS} ${DCA_KERNEL_LIBS}
  )
  #;mc_kernels;lapack;lapack_kernels;blas_kernels;gpu_utils;magma::magma;magma::sparse;

dca_add_gtest(tp_accumulator_multiband_gpu_test
  FAST
  CUDA
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     ${DCA_LIBS} ${DCA_KERNEL_LIBS}
  )

add_subdirectory(ndft)
