//-*-C++-*-

#include <complex>
#include <iostream>
#include <stdexcept>
#include <cstdio>

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

#include "src/linalg_operations/memory_management_GPU.cu.h"

#include "src/linalg_operations/copy_from_CPU_GPU.cu.h"
#include "src/linalg_operations/copy_from_GPU_CPU.cu.h"
#include "src/linalg_operations/copy_from_GPU_GPU.cu.h"

// CUBLAS 1
#include "src/linalg_operations/BLAS_1_SCALE_GPU.cu.h"
#include "src/linalg_operations/BLAS_1_AXPY_GPU.cu.h"
#include "src/linalg_operations/BLAS_1_COPY_GPU.cu.h"
#include "src/linalg_operations/BLAS_1_SWAP_GPU.cu.h"

// CUBLAS 3
#include "src/linalg_operations/BLAS_3_TRSM_GPU.cu.h"
#include "src/linalg_operations/BLAS_3_GEMM_GPU.cu.h"

// own kernels
#include "src/linalg_operations/GEMD_GPU.cu.h"
#include "src/linalg_operations/BENNET_GPU.cu.h"
#include "src/linalg_operations/LU_MATRIX_OPERATIONS_GPU.cu.h"

// magma
#include "src/linalg_operations/LASET_GPU.cu.h"
#include "src/linalg_operations/TRSV_GPU.cu.h"
#include "src/linalg_operations/GETRF_GPU.cu.h"
#include "src/linalg_operations/GETRI_GPU.cu.h"
#include "src/linalg_operations/GETRS_GPU.cu.h"
#include "src/linalg_operations/GESVD_GPU.cu.h"
#include "src/linalg_operations/GEEV_GPU.cu.h"
