// Collects all GPU code.

// #define NDEBUG  // This is now added by CMake to the CUDA_NVCC_FLAGS.
// #define DEBUG_CUDA
// #define cudaDeviceScheduleBlockingSync 0x04

#include "comp_library/linalg/include_linalg.cu.h"

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker/ctaux_walker_routines_GPU.cu.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker/ctaux_walker_tools/ctaux_G0_matrix_routines/ctaux_G0_matrix_routines_GPU.cu.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker/ctaux_walker_tools/ctaux_G_matrix_routines/ctaux_G_matrix_routines_GPU.cu.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker/ctaux_walker_tools/ctaux_N_matrix_routines/ctaux_N_matrix_routines_GPU.cu.h"
