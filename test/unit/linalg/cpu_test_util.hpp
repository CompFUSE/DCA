// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides some utilities to test simple Matrix<CPU> operations.

#ifndef DCA_TEST_UNIT_LINALG_CPU_TEST_UTIL_HPP
#define DCA_TEST_UNIT_LINALG_CPU_TEST_UTIL_HPP

#include "dca/linalg/matrix.hpp"

namespace testing {
// The elements of the matrix will be set with mat(i, j) = func(i, j).
// In: func
// Out: mat
template <typename ScalarType, typename F>
void setMatrixElements(dca::linalg::Matrix<ScalarType, dca::linalg::CPU>& mat, F& func) {
  for (int j = 0; j < mat.nrCols(); ++j)
    for (int i = 0; i < mat.nrRows(); ++i) {
      ScalarType el(func(i, j));
      mat(i, j) = el;
    }
}
}  // testing

#endif  // DCA_TEST_UNIT_LINALG_CPU_TEST_UTIL_HPP
