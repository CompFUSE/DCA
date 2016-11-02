// Collects all GPU code.

// #define NDEBUG  // This is now added by CMake to the CUDA_NVCC_FLAGS.
// #define DEBUG_CUDA
// #define cudaDeviceScheduleBlockingSync 0x04

#include "cuda_runtime.h"

#include "dca/util/integer_division.hpp"
#include "dca/linalg/util/error_cuda.hpp"
#include "dca/linalg/util/stream_functions.hpp"

#include "dca/phys/dca_step/cluster_solver/ctaux/walker/ct_aux_walker_tools.cu.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/g0_interpolation/g0_interpolation.cu.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/g_matrix_tools/g_matrix_tools.cu.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/n_matrix_tools/n_matrix_tools.cu.hpp"
