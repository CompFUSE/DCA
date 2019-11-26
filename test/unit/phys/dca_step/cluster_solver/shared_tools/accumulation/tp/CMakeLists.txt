# test/unit/phys/dca_step/cluster_solver/shared_tools/tp

dca_add_gtest(tp_accumulator_test
  FAST
  GTEST_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     ${DCA_LIBS}
  )

dca_add_gtest(tp_accumulator_singleband_test
  FAST
  GTEST_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     ${DCA_LIBS}
  )

dca_add_gtest(tp_accumulator_particle_hole_test
  FAST
  CUDA
  GTEST_MAIN
  INCLUDE_DIRS ${PROJECT_SOURCE_DIR} ${FFTW_INCLUDE_DIRS} ${HDF5_INCLUDE_DIR}
  LIBS function cluster_domains time_and_frequency_domains quantum_domains timer json random
       enumerations dca_algorithms ${DCA_THREADING_LIBS} parallel_util ${LAPACK_LIBRARIES})

dca_add_gtest(tp_accumulator_gpu_test
  CUDA
  GTEST_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     ${DCA_LIBS}
  )

dca_add_gtest(tp_accumulator_singleband_gpu_test
  CUDA
  GTEST_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS};${PROJECT_SOURCE_DIR}
  LIBS     ${DCA_LIBS}
  )

add_subdirectory(ndft)
