# SS-CT-HYB

add_library(ss_ct_hyb STATIC structures/ss_ct_hyb_configuration.cpp)
target_include_directories(ss_ct_hyb PUBLIC ${FFTW_INCLUDE_DIR})
