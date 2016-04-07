################################################################################
# applications/cluster_solver_check/test/CT-AUX/bilayer_lattice/Nc1_intraband/\
# CMakeLists.txt
################################################################################

file(COPY input.bilayer_lattice_Nc1_intraband.json
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY data.ED.hdf5
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY check_data.QMC.hdf5
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(TEST_INCLUDES "${PROJECT_SOURCE_DIR}/applications;${CMAKE_CURRENT_SOURCE_DIR}/..;${DCA_INCLUDES}")

add_gtest(NAME bilayer_lattice_Nc1_intraband
  ADVANCED
  MPI MPI_NUMPROC 8
        INCLUDES "${TEST_INCLUDES}"
  LIBS     "${DCA_LIBRARIES}"
  )
