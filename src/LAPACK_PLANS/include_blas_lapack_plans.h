//-*-C++-*-

/*!
 *  \author Peter Staar
 */

enum    matrix_form  {GENERAL, HERMITIAN, LU, L, U};
typedef matrix_form  matrix_form_type;

enum    LINEAR_ALGEBRA_LIBRARY {BLAS_LIBRARY, CUBLAS_LIBRARY, LAPACK_LIBRARY, MAGMA_LIBRARY};
typedef LINEAR_ALGEBRA_LIBRARY LINEAR_ALGEBRA_LIBRARY_TYPE;

// BLAS
#include "BLAS/BLAS_C_wrappers.h" 

//#include "swap_plan.h"
//#include "copy_plan.h"
//#include "scale_plan.h"

#include "BLAS/blas_gemm.h"

#include "BLAS/gemv_plan.h"
//#include "gbmv_plan.h"
//#include "trsv_plan.h"

//LAPACK
#include "LAPACK/LAPACK_C_wrappers.h" 

//#include "geev_plan.h"
//#include "include_eigensystem_plans.h"
//#include "include_gelss_plan.h"
#include "LAPACK/geinv_plan.h"
//#include "gelss_plan.h"
#include "LAPACK/gesv_plan.h"
//#include "gesvd_plan.h"
//#include "getrf_plan.h"
//#include "lapack_getrs.h"
#include "LAPACK/geqr_plans/geqr_plans_real.h"

// LINEAR_ALGEBRA_PLANS
#include "LINEAR_ALGEBRA_PLANS/gemm_plan.h"
//#include "getrs_plan.h"
