# test/unit/phys/dca_step/cluster_solver
if (DCA_HAVE_GPU)
  set(KERNELS_LIB mc_kernels)
endif()

add_subdirectory(ctaux)
add_subdirectory(ctint)
add_subdirectory(shared_tools)
add_subdirectory(thread_qmci)
