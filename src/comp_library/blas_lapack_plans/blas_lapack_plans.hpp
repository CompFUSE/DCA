// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file pulls in all BLAS and LAPACK plans *.h files.
// It is self-contained and can be included whenever a file depends on any of the BLAS and LAPACK
// plans *.h files.
//
// TODO: Make all header files self-contained.

#ifndef COMP_LIBRARY_BLAS_LAPACK_PLANS_BLAS_LAPACK_PLANS_HPP
#define COMP_LIBRARY_BLAS_LAPACK_PLANS_BLAS_LAPACK_PLANS_HPP

// BLAS
#include "dca/linalg/blas/blas1.hpp"
#include "dca/linalg/blas/blas2.hpp"
#include "dca/linalg/blas/blas3.hpp"

// LAPACK
#include "comp_library/blas_lapack_plans/LAPACK/LAPACK_C_wrappers.h"
#include "comp_library/blas_lapack_plans/LAPACK/gesv_plan.h"

#endif  // COMP_LIBRARY_BLAS_LAPACK_PLANS_BLAS_LAPACK_PLANS_HPP
