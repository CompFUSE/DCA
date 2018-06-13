// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides reference implementations for Matrix<CPU> operations.

#ifndef DCA_TEST_UNIT_LINALG_MATRIXOP_REFERENCE_HPP
#define DCA_TEST_UNIT_LINALG_MATRIXOP_REFERENCE_HPP

#include "dca/linalg/matrixop.hpp"
#include <complex>
#include <stdexcept>
#include <string>
#include <utility>
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/vector.hpp"

namespace testing {
template <typename ScalarType>
ScalarType conjugate(ScalarType x) {
  return x;
}

template <typename ScalarType>
std::complex<ScalarType> conjugate(std::complex<ScalarType> x) {
  return std::conj(x);
}

// Reference implementation of the matrix-vector multiplication
// In/Out: c ('In' only if beta != 0)
// Preconditions: transa should be one of the following: 'N', 'T' or 'C',
//                a.nrRows() == y.size() if transa == 'N', a.nrCols() == y.size() otherwise,
//                a.nrCols() == x.size() if transa == 'N', a.nrRows() == x.size() otherwise.
// Remark: This implementation is inefficient.
//         It should only be used for testing purpose with small matrices.
template <typename ScalarType>
void refGemv(char transa, ScalarType alpha,
             const dca::linalg::Matrix<ScalarType, dca::linalg::CPU>& a,
             const dca::linalg::Vector<ScalarType, dca::linalg::CPU>& x, ScalarType beta,
             dca::linalg::Vector<ScalarType, dca::linalg::CPU>& y) {
  // Set the values for transa equal 'N'.
  int ma = a.nrRows();
  int na = a.nrCols();

  auto op_a = [transa, &a](int i, int j) {
    switch (transa) {
      case 'T':
        return a(j, i);
      case 'C':
        return testing::conjugate(a(j, i));
      default:
        return a(i, j);
    }
  };

  if (transa == 'T' || transa == 'C') {
    ma = a.nrCols();
    na = a.nrRows();
  }
  else if (transa != 'N') {
    // transa is not 'N', 'T' or 'C'
    throw std::logic_error("Wrong value for transa");
  }

  if (ma != y.size() || na != x.size())
    throw std::logic_error("Wrong matrix sizes");

  for (int i = 0; i < y.size(); ++i) {
    y[i] *= beta;
    for (int j = 0; j < x.size(); ++j)
      y[i] += alpha * op_a(i, j) * x[j];
  }
}

// Reference implementation of the matrix-matrix multiplication
// In/Out: c ('In' only if beta != 0)
// Preconditions: transa and transb should be one of the following: 'N', 'T', 'C',
//                a.nrRows() == c.nrRows() if transa == 'N', a.nrCols() == c.nrRows() otherwise,
//                b.nrCols() == c.nrCols() if transb == 'N', b.nrRows() == c.nrCols() otherwise,
//                ka == kb, where ka = a.nrCols() if transa == 'N', ka = a.nrRows() otherwise and
//                          kb = b.nrRows() if transb == 'N', kb = b.nrCols() otherwise.
// Remark: This implementation is inefficient.
//         It should only be used for testing purpose with small matrices.
template <typename ScalarTypeA, typename ScalarTypeB, typename ScalarTypeC>
void refGemm(char transa, char transb, ScalarTypeC alpha,
             const dca::linalg::Matrix<ScalarTypeA, dca::linalg::CPU>& a,
             const dca::linalg::Matrix<ScalarTypeB, dca::linalg::CPU>& b, ScalarTypeC beta,
             dca::linalg::Matrix<ScalarTypeC, dca::linalg::CPU>& c) {
  static_assert(std::is_same<ScalarTypeA, ScalarTypeC>::value ||
                    std::is_same<std::complex<ScalarTypeA>, ScalarTypeC>::value,
                "Wrong ScalarType for a");
  static_assert(std::is_same<ScalarTypeB, ScalarTypeC>::value ||
                    std::is_same<std::complex<ScalarTypeB>, ScalarTypeC>::value,
                "Wrong ScalarType for b");
  // Set the values for transa and transb equal 'N'.
  int ma = a.nrRows();
  int ka = a.nrCols();
  int kb = b.nrRows();
  int nb = b.nrCols();
  auto op_a = [transa, &a](int i, int j) {
    switch (transa) {
      case 'T':
        return a(j, i);
      case 'C':
        return testing::conjugate(a(j, i));
      default:
        return a(i, j);
    }
  };
  auto op_b = [transb, &b](int i, int j) {
    switch (transb) {
      case 'T':
        return b(j, i);
      case 'C':
        return testing::conjugate(b(j, i));
      default:
        return b(i, j);
    }
  };

  if (transa == 'T' || transa == 'C') {
    ma = a.nrCols();
    ka = a.nrRows();
  }
  else if (transa != 'N') {
    // transa is not 'N', 'T' or 'C'
    throw std::logic_error("Wrong value for transa");
  }
  if (transb == 'T' || transb == 'C') {
    kb = b.nrCols();
    nb = b.nrRows();
  }
  else if (transb != 'N') {
    // transb is not 'N', 'T' or 'C'
    throw std::logic_error("Wrong value for transb");
  }

  if (ma != c.nrRows() || ka != kb || nb != c.nrCols())
    throw std::logic_error("Wrong matrix sizes");

  for (int j = 0; j < c.nrCols(); ++j)
    for (int i = 0; i < c.nrRows(); ++i) {
      c(i, j) *= beta;
      for (int k = 0; k < ka; ++k)
        c(i, j) += alpha * op_a(i, k) * op_b(k, j);
    }
}
}  // testing

#endif  // DCA_TEST_UNIT_LINALG_MATRIXOP_REFERENCE_HPP
