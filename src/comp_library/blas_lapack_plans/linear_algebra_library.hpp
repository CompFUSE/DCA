// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file defines the supported Linear Algebra library types.

#ifndef COMP_LIBRARY_BLAS_LAPACK_PLANS_LINEAR_ALGEBRA_LIBRARY_HPP
#define COMP_LIBRARY_BLAS_LAPACK_PLANS_LINEAR_ALGEBRA_LIBRARY_HPP

enum LINEAR_ALGEBRA_LIBRARY { BLAS_LIBRARY, CUBLAS_LIBRARY, LAPACK_LIBRARY, MAGMA_LIBRARY };
typedef LINEAR_ALGEBRA_LIBRARY LINEAR_ALGEBRA_LIBRARY_TYPE;

#endif  // COMP_LIBRARY_BLAS_LAPACK_PLANS_LINEAR_ALGEBRA_LIBRARY_HPP
