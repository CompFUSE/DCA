# statistical_testing library performance test.

add_executable(stat_test_benchmark statisitical_test_benchmark.cpp)
target_link_libraries(stat_test_benchmark PRIVATE statistical_testing random ${LAPACK_LIBRARIES})
target_include_directories(stat_test_benchmark PRIVATE  ${DCA_INCLUDE_DIRS})
