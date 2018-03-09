# Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)

dca_add_library(statistical_testing STATIC statistical_testing.cpp)

target_include_directories(statistical_testing PRIVATE ${DCA_INCLUDES})
target_link_libraries(statistical_testing  PUBLIC function PRIVATE ${LAPACK_LIBRARIES})
