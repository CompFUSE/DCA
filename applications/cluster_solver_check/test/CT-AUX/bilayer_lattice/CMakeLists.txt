################################################################################
# applications/cluster_solver_check/test/CT-AUX/bilayer_lattice/CMakeLists.txt
################################################################################

file(COPY input.bilayer_lattice.json
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY data.ED.hdf5
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

list(APPEND DCA_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}")

add_gtest(NAME bilayer_lattice
  ADVANCED
  MPI MPI_NUMPROC 8
  INCLUDES "${DCA_INCLUDES}"
  LIBS     "${DCA_LIBRARIES}"
  )
