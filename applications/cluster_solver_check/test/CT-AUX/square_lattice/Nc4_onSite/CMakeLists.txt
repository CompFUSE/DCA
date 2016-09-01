################################################################################
# applications/cluster_solver_check/test/CT-AUX/square_lattice/Nc4_onSite/\
# CMakeLists.txt
################################################################################

dca_add_gtest(square_lattice_Nc4_onSite
  EXTENSIVE
  MPI MPI_NUMPROC 8
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS}
  LIBS ${DCA_LIBS})
