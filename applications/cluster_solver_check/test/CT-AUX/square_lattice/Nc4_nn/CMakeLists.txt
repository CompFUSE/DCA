################################################################################
# applications/cluster_solver_check/test/CT-AUX/square_lattice/Nc4_nn/\
# CMakeLists.txt
################################################################################

dca_add_gtest(square_lattice_Nc4_nn
  EXTENSIVE
  MPI MPI_NUMPROC 8
  PTHREADS
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS}
  LIBS ${DCA_LIBS})
