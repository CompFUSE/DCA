# Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)

if (DCA_HAVE_MULDEFS)
  add_library(statistical_testing STATIC statistical_testing.cpp)
else()
  # We need to create a shared library to avoid a multiple definition error when linking.
  add_library(statistical_testing SHARED statistical_testing.cpp)
endif()

target_include_directories(statistical_testing PRIVATE ${DCA_INCLUDES})
target_link_libraries(statistical_testing  PUBLIC function PRIVATE ${LAPACK_LIBRARIES})
