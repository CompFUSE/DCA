# parallel no_concurrency
set(this_target parallel_no_concurrency)
add_library(${this_target} STATIC no_concurrency.cpp)
dca_gpu_runtime_link(${this_target})
if (DCA_HAVE_ADIOS2)
  target_link_libraries(${this_target} PRIVATE dca_io)
endif ()
