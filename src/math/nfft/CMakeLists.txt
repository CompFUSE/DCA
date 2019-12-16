# NFFT

add_library(nfft STATIC
  domains/nfft_time_domain_names.cpp
  window_functions/gaussian_window_function.cpp
  window_functions/kaiser_bessel_function.cpp)

if (DCA_HAVE_CUDA)
  CUDA_ADD_LIBRARY(dnfft_kernels
    dnfft_1d_kernels.cu
    nfft_helper.cu)
  target_link_libraries(dnfft_kernels mc_kernels)
  target_compile_definitions(dnfft_kernels PUBLIC DCA_HAVE_CUDA)
endif()
