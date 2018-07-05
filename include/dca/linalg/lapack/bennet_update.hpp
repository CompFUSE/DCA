// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides Bennet update implementations.

#ifndef DCA_LINALG_LAPACK_BENNET_UPDATE_HPP
#define DCA_LINALG_LAPACK_BENNET_UPDATE_HPP

#include <complex>

namespace dca {
namespace linalg {
namespace lapack {
// dca::linalg::lapack::

// Peter Stange, Andreas Griewank and Matthias Bollho,
// On the efficient update of rectangular LU factorizations subject to low rank modifications.
template <typename ScalarType>
void rowwiseBennet(int n, int ldm, ScalarType* m, ScalarType* c, ScalarType* r) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      c[i] = c[i] - c[j] * m[i + j * ldm];
      m[i + j * ldm] = m[i + j * ldm] + r[j] * c[i];
    }

    m[i + i * ldm] = m[i + i * ldm] + c[i] * r[i];
    r[i] = r[i] / m[i + i * ldm];

    for (int j = i + 1; j < n; j++) {
      m[i + j * ldm] = m[i + j * ldm] + c[i] * r[j];
      r[j] = r[j] - r[i] * m[i + j * ldm];
    }
  }
}

template <typename ScalarType>
static void standardBennet(int n, int ldm, ScalarType* m, ScalarType* c, ScalarType* r) {
  for (int i = 0; i < n; ++i) {
    m[i + i * ldm] += c[i] * r[i];
    r[i] /= m[i + i * ldm];

    for (int j = i + 1; j < n; ++j) {
      c[j] -= c[i] * m[j + i * ldm];
      m[j + i * ldm] += c[j] * r[i];
    }

    for (int j = i + 1; j < n; ++j) {
      m[i + j * ldm] += c[i] * r[j];
      r[j] -= r[i] * m[i + j * ldm];
    }
  }
}

}  // lapack
}  // linalg
}  // dca

#endif  // DCA_LINALG_LAPACK_BENNET_UPDATE_HPP
