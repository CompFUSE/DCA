# Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)

add_library(statistical_testing STATIC statistical_testing.cpp)

target_include_directories(statistical_testing PRIVATE ${DCA_INCLUDE_DIRS})
target_link_libraries(statistical_testing  PUBLIC function PRIVATE ${LAPACK_LIBRARIES} ${DCA_CUDA_LIBS})
