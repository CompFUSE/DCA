# parallel stdthread
add_library(parallel_stdthread STATIC stdthread.cpp
                                      thread_pool/thread_pool.cpp thread_pool/affinity.cpp)
target_link_libraries(parallel_stdthread pthread)
