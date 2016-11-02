# SS-CT-HYB

if (DCA_HAVE_MULDEFS)
add_library(ss_ct_hyb STATIC structures/ss_ct_hyb_configuration.cpp)

else()
  add_library(ss_ct_hyb SHARED structures/ss_ct_hyb_configuration.cpp)
  target_link_libraries(ss_ct_hyb PUBLIC function)
endif()

target_include_directories(ss_ct_hyb PUBLIC ${FFTW_INCLUDE_DIR})
