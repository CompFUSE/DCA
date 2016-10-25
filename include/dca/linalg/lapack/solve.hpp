// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides the implementation of the solution of system of equations.

#ifndef DCA_LINALG_LAPACK_SOLVE_HPP
#define DCA_LINALG_LAPACK_SOLVE_HPP

#include <complex>
#include "dca/linalg/lapack/lapack.hpp"

namespace dca {
namespace linalg {
namespace lapack {
// dca::linalg::lapack::

// Solves the system of linear equations a * x = b
// using the inverse using LU decomposition.
// Precondition: lda >= n, ldb >= n.
// Postcondition: The matrix a is overwritten by its LU factors.
//                The rhs vectors are overwritten with the solutions of the system.
template <typename Type>
inline void solve(int n, int nrhs, Type* a, int lda, Type* b, int ldb) {
  std::vector<int> ipiv(n);
  gesv(n, nrhs, a, lda, &ipiv[0], b, ldb);
}
template <typename Type>
inline void solve(int n, Type* a, int lda, Type* b) {
  solve(n, 1, a, lda, b, n);
}

}  // lapack
}  // linalg
}  // dca

#endif  // DCA_LINALG_LAPACK_SOLVE_HPP
