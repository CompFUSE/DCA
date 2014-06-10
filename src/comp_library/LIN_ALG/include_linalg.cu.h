//-*-C++-*-

#include "linalg_device_types.h"

#include "cuda.h"
#include "cuda_runtime.h"
#include "cuComplex.h"

#include "cublas.h"
#include "cublas_v2.h"

#define HAVE_CUBLAS 1
#include "magma.h"

#include "basic_cuda_functions.h"
#include "basic_cublas_functions.h"

//#include <thrust/host_vector.h>
//#include <thrust/device_vector.h>

#include <utility>

#include "memory_management_GPU.cu.h"

#include "copy_from_CPU_GPU.cu.h"
#include "copy_from_GPU_CPU.cu.h"
#include "copy_from_GPU_GPU.cu.h"

// cublas
#include "DOT_GPU.cu.h"

// CUBLAS 1
#include "BLAS_1_SCALE_GPU.cu.h"
#include "BLAS_1_AXPY_GPU.cu.h"
#include "BLAS_1_COPY_GPU.cu.h"
#include "BLAS_1_SWAP_GPU.cu.h"

// CUBLAS 3
#include "BLAS_3_TRSM_GPU.cu.h"
#include "BLAS_3_GEMM_GPU.cu.h"

// own kernels
#include "GEMD_GPU.cu.h"
#include "BENNET_GPU.cu.h"
#include "LU_MATRIX_OPERATIONS_GPU.cu.h"

// magma
#include "LASET_GPU.cu.h"
#include "TRSV_GPU.cu.h"
#include "GETRF_GPU.cu.h"
#include "GETRI_GPU.cu.h"
#include "GETRS_GPU.cu.h"
#include "GESVD_GPU.cu.h"
#include "GEEV_GPU.cu.h"

