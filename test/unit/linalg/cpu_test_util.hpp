// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides some utilities to test simple Matrix<CPU> operations.

#ifndef DCA_TEST_UNIT_LINALG_CPU_TEST_UTIL_HPP
#define DCA_TEST_UNIT_LINALG_CPU_TEST_UTIL_HPP

#include "dca/linalg/matrix.hpp"
#include "dca/linalg/vector.hpp"

namespace testing {
// The elements of the matrix will be set with mat(i, j) = func(i, j).
// In: func
// Out: mat
template <typename Matrix, typename F>
void setMatrixElements(Matrix& mat, F&& func) {
  for (int j = 0; j < mat.nrCols(); ++j)
    for (int i = 0; i < mat.nrRows(); ++i) {
      mat(i, j) = func(i, j);
    }
}

// The elements of the vector will be set with vec[i] = func(i).
// In: func
// Out: vec
template <typename Vector, typename F>
void setVectorElements(Vector& vec, F& func) {
  for (int i = 0; i < vec.size(); ++i) {
    vec[i] = func(i);
  }
}
}  // testing

#endif  // DCA_TEST_UNIT_LINALG_CPU_TEST_UTIL_HPP
