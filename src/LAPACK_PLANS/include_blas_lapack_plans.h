//-*-C++-*-

/*!
 *  \author Peter Staar
 */

enum    matrix_form  {GENERAL, HERMITIAN, LU, L, U};
typedef matrix_form  matrix_form_type;

enum    LINEAR_ALGEBRA_LIBRARY {BLAS_LIBRARY, CUBLAS_LIBRARY, LAPACK_LIBRARY, MAGMA_LIBRARY};
typedef LINEAR_ALGEBRA_LIBRARY LINEAR_ALGEBRA_LIBRARY_TYPE;

// BLAS
#include "BLAS_C_wrappers.h" 

//#include "swap_plan.h"
//#include "copy_plan.h"
//#include "scale_plan.h"

#include "blas_gemm.h"

#include "gemv_plan.h"
//#include "gbmv_plan.h"
//#include "trsv_plan.h"

//LAPACK
#include "LAPACK_C_wrappers.h" 

//#include "geev_plan.h"
//#include "include_eigensystem_plans.h"
//#include "include_gelss_plan.h"
#include "geinv_plan.h"
//#include "gelss_plan.h"
#include "gesv_plan.h"
//#include "gesvd_plan.h"
//#include "getrf_plan.h"
//#include "lapack_getrs.h"
#include "geqr_plans_real.h"

// LINEAR_ALGEBRA_PLANS
#include "gemm_plan.h"
//#include "getrs_plan.h"
