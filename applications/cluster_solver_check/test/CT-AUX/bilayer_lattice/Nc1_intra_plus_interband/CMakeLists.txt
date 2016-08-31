################################################################################
# applications/cluster_solver_check/test/CT-AUX/bilayer_lattice/\
# Nc1_intra_plus_interband/CMakeLists.txt
################################################################################

add_gtest(NAME bilayer_lattice_Nc1_intra_plus_interband
  ADVANCED
  MPI MPI_NUMPROC 8
  INCLUDES ${DCA_INCLUDES}
  LIBS ${DCA_LIBRARIES})

target_compile_definitions(bilayer_lattice_Nc1_intra_plus_interband PRIVATE DCA_WITH_REDUCED_VERTEX_FUNCTION)
