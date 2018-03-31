# NFFT

add_library(nfft STATIC
  domains/nfft_time_domain_names.cpp
  window_functions/gaussian_window_function.cpp
  window_functions/kaiser_bessel_function.cpp)
