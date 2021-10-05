# parallel no_concurrency
set(this_target parallel_no_concurrency)
add_library(${this_target} STATIC no_concurrency.cpp)
dca_gpu_runtime_link(${this_target})
