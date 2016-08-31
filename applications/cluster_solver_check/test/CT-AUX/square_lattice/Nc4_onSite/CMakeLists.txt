################################################################################
# applications/cluster_solver_check/test/CT-AUX/square_lattice/Nc4_onSite/\
# CMakeLists.txt
################################################################################

add_gtest(NAME square_lattice_Nc4_onSite
  ADVANCED
  MPI MPI_NUMPROC 8
  INCLUDES ${DCA_INCLUDES}
  LIBS ${DCA_LIBRARIES})

target_compile_definitions(square_lattice_Nc4_onSite PRIVATE DCA_WITH_REDUCED_VERTEX_FUNCTION)
