// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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
#include <tuple>

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
template <typename ScalarType,  template <typename, DeviceType> class MatrixType>
int getInverseWorkSize(MatrixType<ScalarType, CPU>& mat) {
  assert(mat.is_square());

  ScalarType tmp;
  lapack::getri(mat.nrRows(), mat.ptr(), mat.leadingDimension(), nullptr, &tmp, -1);
  return lapack::util::getWorkSize(tmp);
}
#ifdef DCA_HAVE_CUDA
template <typename ScalarType,  template <typename, DeviceType> class MatrixType>
int getInverseWorkSize(const MatrixType<ScalarType, GPU>& mat) {
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
  ScalarType tmp2;
  lapack::geev(&jobvl, &jobvr, mat.nrRows(), mat.ptr(), mat.leadingDimension(), nullptr, nullptr,
               ld, nullptr, ld, &tmp, -1, &tmp2);
  return lapack::util::getWorkSize(tmp);
}
// Returns optimal lwork and liwork for the symmetric eigensolver.
// In: mat
template <typename ScalarType>
std::tuple<int, int> getEigensolverSymmetricWorkSize(char jobv, char uplo,
                                                     Matrix<ScalarType, CPU>& mat) {
  assert(mat.is_square());

  ScalarType tmp1;
  int tmp2;
  lapack::syevd(&jobv, &uplo, mat.nrRows(), mat.ptr(), mat.leadingDimension(), nullptr, &tmp1, -1,
                &tmp2, -1);
  return std::make_tuple(lapack::util::getWorkSize(tmp1), tmp2);
}

// Returns optimal lwork and liwork for the Hermitian eigensolver.
// In: mat
template <typename ScalarType>
std::tuple<int, int, int> getEigensolverHermitianWorkSize(char jobv, char uplo,
                                                          Matrix<std::complex<ScalarType>, CPU>& mat) {
  assert(mat.is_square());

  std::complex<ScalarType> tmp1;
  ScalarType tmp2;
  int tmp3;
  lapack::heevd(&jobv, &uplo, mat.nrRows(), mat.ptr(), mat.leadingDimension(), nullptr, &tmp1, -1,
                &tmp2, -1, &tmp3, -1);
  return std::make_tuple(lapack::util::getWorkSize(tmp1), lapack::util::getWorkSize(tmp2), tmp3);
}
}  // util
}  // matrixop
}  // linalg
}  // dca

#endif  // DCA_LINALG_UTIL_UTIL_MATRIXOP_HPP
