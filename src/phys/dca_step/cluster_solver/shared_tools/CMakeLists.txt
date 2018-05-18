# phys/dca_step/cluster_solver/shared_tools

if (DCA_HAVE_CUDA)
  CUDA_ADD_LIBRARY(mc_tools_kernels
    accumulation/tp/ndft/ndft_kernels.cu
    accumulation/tp/tp_accumulator_kernels.cu)
endif()
