// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi(gbalduzz@itp.phys.ethz.ch)
//
// This file provides a lightweight proxy to access blocks of a matrix. The underlying matrix must
// not be destroyed while this object is in use. A const and a non const version are provided.

#ifndef DCA_LINALG_MAKE_CONSTANT_VIEW_HPP
#define DCA_LINALG_MAKE_CONSTANT_VIEW_HPP

#include <memory>

#include "dca/linalg/matrix_view.hpp"

namespace dca {
namespace linalg {
// dca::linalg::

template <typename ScalarType, DeviceType device_t>
auto inline makeConstantView(const ScalarType* data, const std::pair<int, int> size, const int ld) {
  return std::make_unique<const MatrixView<ScalarType, device_t>>(const_cast<ScalarType*>(data),
                                                                  size, ld);
}

template <typename ScalarType, DeviceType device_t>
inline auto makeConstantView(ScalarType* const data, const int size, const int ld) {
  return makeConstantView<ScalarType, device_t>(data, std::make_pair(size, size), ld);
}

template <typename ScalarType, DeviceType device_t>
inline auto makeConstantView(ScalarType* const data, const int size) {
  return makeConstantView<ScalarType, device_t>(data, std::make_pair(size, size), size);
}

template <typename ScalarType, DeviceType device_t, template <typename, DeviceType> class Matrix>
inline auto makeConstantView(const Matrix<ScalarType, device_t>& mat) {
  return makeConstantView<ScalarType, device_t>(mat.ptr(), mat.size(), mat.leadingDimension());
}

template <typename ScalarType, DeviceType device_t, template <typename, DeviceType> class Matrix>
inline auto makeConstantView(const Matrix<ScalarType, device_t>& mat, const int offset_i,
                             const int offset_j) {
  assert(offset_i < mat.nrCols());
  assert(offset_j < mat.nrRows());
  return makeConstantView<ScalarType, device_t>(
      mat.ptr(offset_i, offset_j), std::make_pair(mat.nrRows() - offset_i, mat.nrCols() - offset_j),
      mat.leadingDimension());
}

template <typename ScalarType, DeviceType device_t, template <typename, DeviceType> class Matrix>
inline auto makeConstantView(const Matrix<ScalarType, device_t>& mat, const int offset_i,
                             const int offset_j, const int ni, const int nj) {
  assert(ni + offset_i <= mat.nrRows());
  assert(nj + offset_j <= mat.nrCols());
  return makeConstantView<ScalarType, device_t>(mat.ptr(offset_i, offset_j), std::make_pair(ni, nj),
                                                mat.leadingDimension());
}

}  // linalg
}  // dca

#endif  // DCA_LINALG_MAKE_CONSTANT_VIEW_HPP
