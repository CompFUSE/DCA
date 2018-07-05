// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides the conversion between char and cublas types for diag, side, trans and uplo.

#ifndef DCA_LINALG_BLAS_CUBLAS_CONVERSION_CHAR_TYPES_HPP
#define DCA_LINALG_BLAS_CUBLAS_CONVERSION_CHAR_TYPES_HPP

#include <cublas_v2.h>
#include <stdexcept>

namespace dca {
namespace linalg {
namespace cublas {
// dca::linalg::cublas::

// Returns the corresponding cublasDiagType_t value to diag.
// Preconditions: diag == 'N', or diag == 'U'.
inline cublasDiagType_t getCublasDiagValue(char diag) {
  switch (diag) {
    case 'U':
      return CUBLAS_DIAG_UNIT;
    case 'N':
      return CUBLAS_DIAG_NON_UNIT;
  }

  throw std::logic_error(__FUNCTION__);
  return cublasDiagType_t();
}

// Returns the corresponding cublasSideMode_t value to side.
// Preconditions: side == 'L', or side == 'R'.
inline cublasSideMode_t getCublasSideValue(char side) {
  switch (side) {
    case 'L':
      return CUBLAS_SIDE_LEFT;
    case 'R':
      return CUBLAS_SIDE_RIGHT;
  }

  throw std::logic_error(__FUNCTION__);
  return cublasSideMode_t();
}

// Returns the corresponding cublasOperation_t value to trans.
// Preconditions: trans == 'N', trans == 'T', or trans == 'C'.
inline cublasOperation_t getCublasTransValue(char trans) {
  switch (trans) {
    case 'N':
      return CUBLAS_OP_N;
    case 'T':
      return CUBLAS_OP_T;
    case 'C':
      return CUBLAS_OP_C;
  }

  throw std::logic_error(__FUNCTION__);
  return cublasOperation_t();
}

// Returns the corresponding cublasFillMode_t value to uplo.
// Preconditions: uplo == 'L', or uplo == 'U'.
inline cublasFillMode_t getCublasUploValue(char uplo) {
  switch (uplo) {
    case 'L':
      return CUBLAS_FILL_MODE_LOWER;
    case 'U':
      return CUBLAS_FILL_MODE_UPPER;
  }

  throw std::logic_error(__FUNCTION__);
  return cublasFillMode_t();
}

}  // cublas
}  // linalg
}  // dca

#endif  // DCA_LINALG_BLAS_CUBLAS_CONVERSION_CHAR_TYPES_HPP
