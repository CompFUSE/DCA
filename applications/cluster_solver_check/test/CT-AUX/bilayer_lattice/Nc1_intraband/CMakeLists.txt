################################################################################
# applications/cluster_solver_check/test/CT-AUX/bilayer_lattice/Nc1_intraband/\
# CMakeLists.txt
################################################################################

add_gtest(NAME bilayer_lattice_Nc1_intraband
          ADVANCED
          MPI MPI_NUMPROC 8
          INCLUDES ${DCA_INCLUDES}
          LIBS ${DCA_LIBRARIES})
