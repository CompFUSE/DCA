# phys/dca_step/cluster_solver/shared_tools

if (DCA_HAVE_CUDA)
  CUDA_ADD_LIBRARY(mc_tools_kernels
    accumulation/tp/ndft/ndft_kernels.cu
    accumulation/tp/tp_accumulator_kernels.cu
    accumulation/tp/g4_helper.cu)
  set_target_properties(mc_tools_kernels PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#  set_property(TARGET mc_tools_kernels PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()
