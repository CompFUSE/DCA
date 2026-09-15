add_library(main_parameters STATIC main_parameters.cpp)
target_include_directories(main_parameters PRIVATE ${DCA_INCLUDE_DIRS} ${SIMPLEX_GM_RULE_INCLUDE_DIR})
# somehow Parameters is dependent on locating fftw3.h!
target_link_libraries(main_parameters PRIVATE FFTW::Double dca_io quantum_domains)
