################################################################################
# applications/cluster_solver_check/test/CT-AUX/bilayer_lattice/Nc1_intraband/\
# CMakeLists.txt
################################################################################

dca_add_gtest(bilayer_lattice_Nc1_intraband
  EXTENSIVE
  MPI MPI_NUMPROC 8
  INCLUDE_DIRS ${DCA_INCLUDES}
  LIBS ${DCA_LIBRARIES})
