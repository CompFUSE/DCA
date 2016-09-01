################################################################################
# applications/cluster_solver_check/test/CT-AUX/square_lattice/\
# Nc4_onSite_plus_nn/CMakeLists.txt
################################################################################

dca_add_gtest(square_lattice_Nc4_onSite_plus_nn
  EXTENSIVE
  MPI MPI_NUMPROC 8
  INCLUDE_DIRS ${DCA_INCLUDES}
  LIBS ${DCA_LIBRARIES})
