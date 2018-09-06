dca_add_gtest(square_lattice_Nc4_onSite_plus_nn
  EXTENSIVE
  MPI MPI_NUMPROC 8
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS}
  LIBS ${DCA_LIBS} ${DCA_CUDA_LIBS})
