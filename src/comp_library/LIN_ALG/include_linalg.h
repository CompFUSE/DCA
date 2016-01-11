//-*-C++-*-

// BLAS
//#include "C_wrappers_BLAS.h" 

//LAPACK
//#include "C_wrappers_LAPACK.h" 

#include "linalg_device_types.h"

#include "src/matrix_scalartype.h"

#include "src/linalg_structures/cublas_thread_manager_tem.h"
#include "src/linalg_structures/cublas_thread_manager_CPU.h"
#include "src/linalg_structures/cublas_thread_manager_GPU.h"

#include "src/linalg_operations/copy_from_tem.h"
#include "src/linalg_operations/copy_from_CPU_CPU.h"
#include "src/linalg_operations/copy_from_CPU_GPU.h"
#include "src/linalg_operations/copy_from_GPU_CPU.h"
#include "src/linalg_operations/copy_from_GPU_GPU.h"

#include "src/linalg_operations/memory_management_tem.h"
#include "src/linalg_operations/memory_management_CPU.h"
#include "src/linalg_operations/memory_management_GPU.h"

#include "src/vector.h"
#include "src/matrix.h"

#include "src/linalg_operations/LU_MATRIX_OPERATIONS.h"
#include "src/linalg_operations/LU_MATRIX_OPERATIONS_CPU.h"
#include "src/linalg_operations/LU_MATRIX_OPERATIONS_GPU.h"

#include "src/linalg_operations/REMOVE_tem.h"
#include "src/linalg_operations/REMOVE_CPU.h"
#include "src/linalg_operations/REMOVE_GPU.h"

// BLAS 1

#include "src/linalg_operations/BLAS_1_AXPY_tem.h"
#include "src/linalg_operations/BLAS_1_AXPY_CPU.h"
#include "src/linalg_operations/BLAS_1_AXPY_GPU.h"

#include "src/linalg_operations/BLAS_1_COPY_tem.h"
#include "src/linalg_operations/BLAS_1_COPY_CPU.h"
#include "src/linalg_operations/BLAS_1_COPY_GPU.h"

#include "src/linalg_operations/BLAS_1_SCALE_tem.h"
#include "src/linalg_operations/BLAS_1_SCALE_CPU.h"
#include "src/linalg_operations/BLAS_1_SCALE_GPU.h"

#include "src/linalg_operations/BLAS_1_SWAP_tem.h"
#include "src/linalg_operations/BLAS_1_SWAP_CPU.h"
#include "src/linalg_operations/BLAS_1_SWAP_GPU.h"

// BLAS 2

#include "src/linalg_operations/BLAS_2_GEMV_tem.h"
#include "src/linalg_operations/BLAS_2_GEMV_CPU.h"

// BLAS 3

#include "src/linalg_operations/BLAS_3_TRSM_tem.h"
#include "src/linalg_operations/BLAS_3_TRSM_CPU.h"
#include "src/linalg_operations/BLAS_3_TRSM_GPU.h"

#include "src/linalg_operations/BLAS_3_GEMM_tem.h"
#include "src/linalg_operations/BLAS_3_GEMM_CPU.h"
#include "src/linalg_operations/BLAS_3_GEMM_GPU.h"


#include "src/linalg_operations/LASET_tem.h"
#include "src/linalg_operations/LASET_CPU.h"
#include "src/linalg_operations/LASET_GPU.h"

#include "src/linalg_operations/DOT_tem.h"
#include "src/linalg_operations/DOT_CPU.h"
#include "src/linalg_operations/DOT_GPU.h"

#include "src/linalg_operations/GEMD_tem.h"
#include "src/linalg_operations/GEMD_CPU.h"
#include "src/linalg_operations/GEMD_GPU.h"

#include "src/linalg_operations/TRSV_tem.h"
#include "src/linalg_operations/TRSV_CPU.h"
#include "src/linalg_operations/TRSV_GPU.h"


#include "src/linalg_operations/BENNET_tem.h"
#include "src/linalg_operations/BENNET_CPU.h"
#include "src/linalg_operations/BENNET_GPU.h"

#include "src/linalg_operations/GETRS_tem.h"
#include "src/linalg_operations/GETRS_CPU.h"
#include "src/linalg_operations/GETRS_GPU.h"

#include "src/linalg_operations/GETRF_tem.h"
#include "src/linalg_operations/GETRF_CPU.h"
#include "src/linalg_operations/GETRF_GPU.h"

#include "src/linalg_operations/GETRI_tem.h"
#include "src/linalg_operations/GETRI_CPU.h"
#include "src/linalg_operations/GETRI_GPU.h"

#include "src/linalg_operations/GEINV_tem.h"

#include "src/linalg_operations/GEEV_tem.h"
#include "src/linalg_operations/GEEV_CPU.h"
#include "src/linalg_operations/GEEV_GPU.h"

#include "src/linalg_operations/GESV_tem.h"
#include "src/linalg_operations/GESV_CPU.h"

#include "src/linalg_operations/GESVD_tem.h"
#include "src/linalg_operations/GESVD_CPU.h"
#include "src/linalg_operations/GESVD_GPU.h"

#include "src/linalg_operations/PSEUDO_INVERSE_tem.h"
#include "src/linalg_operations/PSEUDO_INVERSE_CPU.h"

// performance_inspector

#include "src/linalg_structures/performance_inspector_tem.h"
#include "src/linalg_structures/performance_inspector_CPU.h"
#include "src/linalg_structures/performance_inspector_GPU.h"
