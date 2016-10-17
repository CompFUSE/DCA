// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file provides the utilities to get workspace size for the matrix operations:
// - inverse,

#ifndef DCA_LINALG_UTIL_UTIL_MATRIXOP_HPP
#define DCA_LINALG_UTIL_UTIL_MATRIXOP_HPP

#include <cassert>
#include <stdexcept>

#include "dca/linalg/matrix.hpp"
#include "dca/linalg/util/util_lapack.hpp"

#ifdef DCA_HAVE_CUDA
#include "dca/linalg/lapack/magma.hpp"
#endif

namespace dca {
namespace linalg {
namespace matrixop {
namespace util {
// dca::linalg::matrixop::

// Returns optimal lwork for inverse.
// In: mat
template <typename ScalarType>
int getInverseWorkSize(Matrix<ScalarType, CPU>& mat) {
  assert(mat.is_square());

  ScalarType tmp;
  lapack::getri(mat.nrRows(), mat.ptr(), mat.leadingDimension(), nullptr, &tmp, -1);
  return lapack::util::getWorkSize(tmp);
}
#ifdef DCA_HAVE_CUDA
template <typename ScalarType>
int getInverseWorkSize(const Matrix<ScalarType, GPU>& mat) {
  assert(mat.is_square());

  return mat.nrRows() * magma::get_getri_nb<ScalarType>(mat.nrRows());
}
#endif  // DCA_HAVE_CUDA

// Returns optimal lwork for the eigensolver.
// In: mat
template <typename ScalarType>
int getEigensolverWorkSize(char jobvl, char jobvr, Matrix<ScalarType, CPU>& mat) {
  assert(mat.is_square());

  int ld = mat.nrRows();
  ScalarType tmp;
  lapack::geev(&jobvl, &jobvr, mat.nrRows(), mat.ptr(), mat.leadingDimension(), nullptr, nullptr,
               nullptr, ld, nullptr, ld, &tmp, -1);
  return lapack::util::getWorkSize(tmp);
}
template <typename ScalarType>
int getEigensolverWorkSize(char jobvl, char jobvr, Matrix<std::complex<ScalarType>, CPU>& mat) {
  assert(mat.is_square());

  int ld = mat.nrRows();
  std::complex<ScalarType> tmp;
  lapack::geev(&jobvl, &jobvr, mat.nrRows(), mat.ptr(), mat.leadingDimension(), nullptr, nullptr,
               ld, nullptr, ld, &tmp, -1, nullptr);
  return lapack::util::getWorkSize(tmp);
}

}  // util
}  // matrixop
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_UTIL_MATRIXOP_HPP
