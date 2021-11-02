# Cluster solver

add_subdirectory(ctaux)
add_subdirectory(ctint)
add_subdirectory(ss_ct_hyb)

if (DCA_HAVE_GPU)
  set(kernel_sources
    # CT-INT
    ctint/device_helper/ctint_helper.cu
    ctint/walker/walker_kernels.cu
    ctint/walker/tools/d_matrix_kernels.cu
    ctint/walker/tools/g0_interpolation_gpu.cu
    # CT-AUX
    ctaux/walker/ct_aux_walker_tools_kernels.cu
    ctaux/walker/tools/g0_interpolation/g0_interpolation_kernels.cu
    ctaux/walker/tools/g_matrix_tools/g_matrix_tools_kernels.cu
    ctaux/walker/tools/n_matrix_tools/n_matrix_tools_kernels.cu
    # Shared Tools
    # Tp tools.
    shared_tools/accumulation/tp/ndft/ndft_kernels.cu
    shared_tools/accumulation/tp/tp_accumulator_kernels.cu
    shared_tools/accumulation/tp/g4_helper.cu
    shared_tools/cluster_helper.cu
    # Sp tools
    ../../../math/nfft/dnfft_1d_kernels.cu
    ../../../math/nfft/nfft_helper.cu
    )
  add_library(mc_kernels ${kernel_sources})

  dca_gpu_runtime_link(mc_kernels)
  dca_gpu_blas_link(mc_kernels)
  dca_gpu_device_link(mc_kernels)
  
  target_link_libraries(mc_kernels PUBLIC ctint)
endif()
