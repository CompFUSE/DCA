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
// TODO: - Make all header files self-contained.
//       - Remove all deprecated files.

#ifndef COMP_LIBRARY_BLAS_LAPACK_PLANS_BLAS_LAPACK_PLANS_HPP
#define COMP_LIBRARY_BLAS_LAPACK_PLANS_BLAS_LAPACK_PLANS_HPP

#include <cassert>
#include <complex>
#include <cstring>
#include <iostream>

enum matrix_form { GENERAL, HERMITIAN, LU, L, U };
using matrix_form_type = matrix_form;

enum LINEAR_ALGEBRA_LIBRARY { BLAS_LIBRARY, CUBLAS_LIBRARY, LAPACK_LIBRARY, MAGMA_LIBRARY };
using LINEAR_ALGEBRA_LIBRARY_TYPE = LINEAR_ALGEBRA_LIBRARY;

// BLAS
#include "BLAS/BLAS_C_wrappers.h"
#include "BLAS/blas_gemm.h"
#include "BLAS/gemv_plan.h"

// LAPACK
#include "LAPACK/LAPACK_C_wrappers.h"
#include "LAPACK/geinv_plan.h"
#include "LAPACK/geqr_plans/geqr_plans_real.h"
#include "LAPACK/gesv_plan.h"

// LINEAR_ALGEBRA_PLANS
#include "LINEAR_ALGEBRA_PLANS/gemm_plan.h"

#endif  // COMP_LIBRARY_BLAS_LAPACK_PLANS_BLAS_LAPACK_PLANS_HPP
