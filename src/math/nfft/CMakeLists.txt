# NFFT

# Need to prefix with 'dca_' to avoid name collision with external NFFT3 library.
add_library(dca_nfft
  STATIC
  domains/nfft_time_domain_names.cpp
  window_functions/gaussian_window_function.cpp)
