dca_add_gtest(stdthread_ctaux_tp_test
  EXTENSIVE
  GTEST_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS}
  LIBS         ${DCA_LIBS}
  )

# nonshared threading not currently support for ctint
# dca_add_gtest(stdthread_ctint_cluster_solver_nonshared_test
#   EXTENSIVE
#   GTEST_MAIN
#   INCLUDE_DIRS ${DCA_INCLUDE_DIRS}
#   LIBS         ${DCA_LIBS}
#   )

dca_add_gtest(stdthread_ctint_cluster_solver_shared_test
  EXTENSIVE
  GTEST_MAIN
  INCLUDE_DIRS ${DCA_INCLUDE_DIRS}
  LIBS         ${DCA_LIBS}
  )

if (DCA_HAVE_ADIOS2)
  if (DCA_HAVE_MPI)
    dca_add_gtest(stdthread_ctint_cluster_solver_io_test
      EXTENSIVE
      INCLUDE_DIRS ${DCA_INCLUDE_DIRS}
      LIBS         ${DCA_LIBS}
      )
  endif()
endif()
add_subdirectory(gpu)
