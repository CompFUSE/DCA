// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// This file defines various matrix forms that are used in the BLAS and LAPACK plans.

#ifndef COMP_LIBRARY_BLAS_LAPACK_PLANS_MATRIX_FORM_HPP
#define COMP_LIBRARY_BLAS_LAPACK_PLANS_MATRIX_FORM_HPP

enum matrix_form { GENERAL, HERMITIAN, LU, L, U };
typedef matrix_form matrix_form_type;

#endif  // COMP_LIBRARY_BLAS_LAPACK_PLANS_MATRIX_FORM_HPP
