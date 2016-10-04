// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file pulls in all Linalg *.h files.
// It is self-contained and can be included whenever a file depends on any of the Linalg *.h files.
//
// TODO: Make all header files self-contained.

#ifndef COMP_LIBRARY_LINALG_LINALG_HPP
#define COMP_LIBRARY_LINALG_LINALG_HPP

#include <cassert>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/linalg/linalg.hpp"

#include "comp_library/blas_lapack_plans/blas_lapack_plans.hpp"

#include "linalg_device_types.h"

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

#include "src/linalg_operations/LU_MATRIX_OPERATIONS.h"
#include "src/linalg_operations/LU_MATRIX_OPERATIONS_CPU.h"
#include "src/linalg_operations/LU_MATRIX_OPERATIONS_GPU.h"

// LAPACK
#include "src/linalg_operations/LASET_tem.h"
#include "src/linalg_operations/LASET_CPU.h"
#include "src/linalg_operations/LASET_GPU.h"

#include "src/linalg_operations/GEMD_tem.h"
#include "src/linalg_operations/GEMD_CPU.h"
#include "src/linalg_operations/GEMD_GPU.h"

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

// Performance inspector
#include "src/linalg_structures/performance_inspector_tem.h"
#include "src/linalg_structures/performance_inspector_CPU.h"
#ifdef DCA_HAVE_CUDA
#include "src/linalg_structures/performance_inspector_GPU.h"
#endif

#endif  // COMP_LIBRARY_LINALG_LINALG_HPP
