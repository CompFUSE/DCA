################################################################################
# applications/cluster_solver_check/test/CT-AUX/bilayer_lattice/Nc1_interband/\
# CMakeLists.txt
################################################################################
message(tets include ${TEST_INCLUDES})

file(COPY input.bilayer_lattice_Nc1_interband.json
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY data.ED.hdf5
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY check_data.QMC.hdf5
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_gtest(NAME bilayer_lattice_Nc1_interband
  ADVANCED
  MPI MPI_NUMPROC 8
  INCLUDES "${TEST_INCLUDES}"
  LIBS     "${DCA_LIBRARIES}"
  )
