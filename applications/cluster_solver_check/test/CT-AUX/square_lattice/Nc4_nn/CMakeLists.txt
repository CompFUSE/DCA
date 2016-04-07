################################################################################
# applications/cluster_solver_check/test/CT-AUX/square_lattice/Nc4_nn/\
# CMakeLists.txt
################################################################################

file(COPY input.square_lattice_Nc4_nn.json
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY data.ED.hdf5
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY check_data.QMC.hdf5
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(TEST_INCLUDES "${PROJECT_SOURCE_DIR}/applications;${CMAKE_CURRENT_SOURCE_DIR}/..;${DCA_INCLUDES}")

add_gtest(NAME square_lattice_Nc4_nn
  ADVANCED
  MPI MPI_NUMPROC 8
  INCLUDES "${TEST_INCLUDES}"
  LIBS     "${DCA_LIBRARIES}"
  )
