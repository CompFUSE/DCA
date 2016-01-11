// Collects all GPU code.

// #define NDEBUG  // This is now added by CMake to the CUDA_NVCC_FLAGS.
// #define DEBUG_CUDA
// #define cudaDeviceScheduleBlockingSync 0x04

#include "include_linalg.cu.h"

#include "ctaux_G0_matrix_routines_GPU.cu.h"
#include "ctaux_G_matrix_routines_GPU.cu.h"
#include "ctaux_N_matrix_routines_GPU.cu.h"

#include "ctaux_walker_routines_GPU.cu.h"
