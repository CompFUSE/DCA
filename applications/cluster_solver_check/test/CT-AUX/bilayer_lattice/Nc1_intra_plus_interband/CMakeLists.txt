################################################################################
# applications/cluster_solver_check/test/CT-AUX/bilayer_lattice/\
# Nc1_intra_plus_interband/CMakeLists.txt
################################################################################

dca_add_gtest(bilayer_lattice_Nc1_intra_plus_interband
  EXTENSIVE
  MPI MPI_NUMPROC 8
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS}
  LIBS ${DCA_LIBS})
