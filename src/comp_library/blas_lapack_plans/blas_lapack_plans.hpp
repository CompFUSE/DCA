// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// This file pulls in all BLAS and LAPACK plans *.h files.
// It is self-contained and can be included whenever a file depends on any of the BLAS and LAPACK
// plans *.h files.
//
// TODO: Make all header files self-contained.

#ifndef COMP_LIBRARY_BLAS_LAPACK_PLANS_BLAS_LAPACK_PLANS_HPP
#define COMP_LIBRARY_BLAS_LAPACK_PLANS_BLAS_LAPACK_PLANS_HPP

// BLAS
#include "comp_library/blas_lapack_plans/BLAS/BLAS_C_wrappers.h"
#include "comp_library/blas_lapack_plans/BLAS/blas_gemm.h"
// #include "comp_library/blas_lapack_plans/BLAS/gbmv_plan.h"
#include "comp_library/blas_lapack_plans/BLAS/gemv_plan.h"

// LAPACK
#include "comp_library/blas_lapack_plans/LAPACK/LAPACK_C_wrappers.h"
// #include "comp_library/blas_lapack_plans/LAPACK/eigensystem_plans/include_eigensystem_plans.h"
#include "comp_library/blas_lapack_plans/LAPACK/geinv_plan.h"
// #include "comp_library/blas_lapack_plans/LAPACK/gelss_plans/include_gelss_plan.h"
#include "comp_library/blas_lapack_plans/LAPACK/geqr_plans/geqr_plans_real.h"
#include "comp_library/blas_lapack_plans/LAPACK/gesv_plan.h"
// #include "comp_library/blas_lapack_plans/LAPACK/gesvd_plan.h"
// #include "comp_library/blas_lapack_plans/LAPACK/getrf_plan.h"
// #include "comp_library/blas_lapack_plans/LAPACK/lapack_getrs.h"

// LINEAR_ALGEBRA_PLANS
#include "comp_library/blas_lapack_plans/LINEAR_ALGEBRA_PLANS/gemm_plan.h"
// #include "comp_library/blas_lapack_plans/LINEAR_ALGEBRA_PLANS/getrs_plan.h"

#endif  // COMP_LIBRARY_BLAS_LAPACK_PLANS_BLAS_LAPACK_PLANS_HPP
