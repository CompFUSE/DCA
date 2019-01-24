dca_add_gtest(bilayer_lattice_Nc1_interband
  EXTENSIVE
  MPI MPI_NUMPROC 8
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS}
  LIBS ${DCA_LIBS} ${DCA_CUDA_LIBS})
