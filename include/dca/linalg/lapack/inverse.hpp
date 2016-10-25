// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides the implementation of the matrix inverse function.

#ifndef DCA_LINALG_LAPACK_INVERSE_HPP
#define DCA_LINALG_LAPACK_INVERSE_HPP

#include <complex>
#include "dca/linalg/lapack/lapack.hpp"

namespace dca {
namespace linalg {
namespace lapack {
// dca::linalg::lapack::

// Computes the inverse using LU decomposition.
// Postcondition: The matrix a is overwritten by its inverse.
template <typename Type>
inline void inverse(int n, Type* a, int lda, int* ipiv, Type* work, int lwork) {
  getrf(n, n, a, lda, ipiv);
  getri(n, a, lda, ipiv, work, lwork);
}
template <typename Type>
inline void inverse(int n, Type* a, int lda) {
  std::vector<int> ipiv(n);
  Type tmp;

  // Get optimal lwork.
  getri(n, a, lda, &ipiv[0], &tmp, -1);
  std::vector<Type> work(util::getWorkSize(tmp));
  inverse(n, a, lda, &ipiv[0], &work[0], work.size());
}

}  // lapack
}  // linalg
}  // dca

#endif  // DCA_LINALG_LAPACK_INVERSE_HPP
