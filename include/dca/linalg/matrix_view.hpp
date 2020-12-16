// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi(gbalduzz@itp.phys.ethz.ch)
//
// This file provides a lightweight proxy to access blocks of a matrix. The underlying matrix must
// not be destroyed while this object is in use. A const and a non const version are provided.

#ifndef DCA_LINALG_MATRIX_VIEW_HPP
#define DCA_LINALG_MATRIX_VIEW_HPP

#include <cassert>
#include <iostream>
#include <memory>
#include <stdexcept>

#include "dca/linalg/device_type.hpp"
#include "dca/util/cuda_definitions.hpp"

namespace dca {
namespace linalg {
// dca::linalg::

template <class Scalar, DeviceType device_name = linalg::CPU>
class MatrixView {
public:
  using ValueType = Scalar;
  constexpr static DeviceType device = device_name;

  using Pair = std::pair<int, int>;
  MatrixView(Scalar* data, Pair size);
  MatrixView(Scalar* data, Pair size, int ld);
  MatrixView(Scalar* data, int size, int ld);
  MatrixView(Scalar* data, int size);

  template <template <class, DeviceType> class Matrix>
  MatrixView(Matrix<Scalar, device_name>& mat);
  template <template <class, DeviceType> class Matrix>
  MatrixView(const Matrix<std::remove_cv_t<Scalar>, device_name>& mat);
  template <template <class, DeviceType> class Matrix>
  MatrixView(Matrix<Scalar, device_name>& mat, int offset_i, int offset_j);
  template <template <class, DeviceType> class Matrix>
  MatrixView(const Matrix<std::remove_cv_t<Scalar>, device_name>& mat, int offset_i, int offset_j);
  template <template <class, DeviceType> class Matrix>
  MatrixView(Matrix<Scalar, device_name>& mat, int offset_i, int offset_j, int ni, int nj);
  template <template <class, DeviceType> class Matrix>
  MatrixView(const Matrix<std::remove_cv_t<Scalar>, device_name>& mat, int offset_i, int offset_j,
             int ni, int nj);

  // Preconditions: Device is CPU.
  //                Assignment from Scalar2 to Scalar is defined.
  //                Size must be equal to rhs' size.
  template <template <class, DeviceType> class Matrix, class Scalar2>
  MatrixView& operator=(const Matrix<Scalar2, device_name>& rhs);

  // Same as above, necessary to override default operator=.
  MatrixView& operator=(const MatrixView& rhs);

  void print(std::ostream& out = std::cout) const;

  __DEVICE__ __HOST__ inline int leadingDimension() const {
    return ldm_;
  }
  __DEVICE__ __HOST__ inline int ld() const {
    return leadingDimension();
  }

  __DEVICE__ __HOST__ std::pair<int, int> size() const {
    return size_;
  }
  __DEVICE__ __HOST__ Scalar* ptr() {
    return ptr_;
  }
  __DEVICE__ __HOST__ const Scalar* ptr() const {
    return ptr_;
  }
  __DEVICE__ __HOST__ Scalar* ptr(int i, int j) {
    assert(0 <= i && i <= size_.first);
    assert(0 <= j && j <= size_.second);
    return ptr_ + leadingDimension() * j + i;
  }
  __DEVICE__ __HOST__ const Scalar* ptr(int i, int j) const {
    assert(0 <= i && i <= size_.first);
    assert(0 <= j && j <= size_.second);
    return ptr_ + leadingDimension() * j + i;
  }
  __DEVICE__ __HOST__ int nrRows() const {
    return size_.first;
  }
  __DEVICE__ __HOST__ int nrCols() const {
    return size_.second;
  }
  bool is_square() const {
    return size_.first == size_.second;
  }
  __DEVICE__ __HOST__ Scalar& operator()(int i, int j) {
    assert(0 <= i && i < size_.first);
    assert(0 <= j && j < size_.second);
    return ptr_[i + j * ldm_];
  }
  __DEVICE__ __HOST__ const Scalar& operator()(int i, int j) const {
    assert(0 <= i && i < size_.first);
    assert(0 <= j && j < size_.second);
    return ptr_[i + j * ldm_];
  }

private:
  Scalar* const ptr_;
  const int ldm_;
  const std::pair<int, int> size_;
};

template <class Scalar, DeviceType device_t>
MatrixView<Scalar, device_t>::MatrixView(Scalar* const data, const Pair size)
    : MatrixView(data, size, size.first) {}

template <class Scalar, DeviceType device_t>
MatrixView<Scalar, device_t>::MatrixView(Scalar* const data, const Pair size, const int ld)
    : ptr_(data), ldm_(ld), size_(size) {}

template <class Scalar, DeviceType device_t>
MatrixView<Scalar, device_t>::MatrixView(Scalar* const data, const int size, const int ld)
    : ptr_(data), ldm_(ld), size_(std::make_pair(size, size)) {}

template <class Scalar, DeviceType device_t>
MatrixView<Scalar, device_t>::MatrixView(Scalar* const data, const int size)
    : ptr_(data), ldm_(size), size_(std::make_pair(size, size)) {}

template <class Scalar, DeviceType device_t>
template <template <class, DeviceType> class Matrix>
MatrixView<Scalar, device_t>::MatrixView(Matrix<Scalar, device_t>& mat)
    : ptr_(mat.ptr()), ldm_(mat.leadingDimension()), size_(mat.size()) {}

template <class Scalar, DeviceType device_t>
template <template <class, DeviceType> class Matrix>
MatrixView<Scalar, device_t>::MatrixView(const Matrix<std::remove_cv_t<Scalar>, device_t>& mat)
    : ptr_(mat.ptr()), ldm_(mat.leadingDimension()), size_(mat.size()) {}

template <class Scalar, DeviceType device_t>
template <template <class, DeviceType> class Matrix>
MatrixView<Scalar, device_t>::MatrixView(Matrix<Scalar, device_t>& mat, int offset_i, int offset_j)
    : MatrixView(mat, offset_i, offset_j, mat.nrRows() - offset_i, mat.nrCols() - offset_j) {
  assert(offset_i < mat.nrCols());
  assert(offset_j < mat.nrRows());
}
template <class Scalar, DeviceType device_t>
template <template <class, DeviceType> class Matrix>
MatrixView<Scalar, device_t>::MatrixView(const Matrix<std::remove_cv_t<Scalar>, device_t>& mat,
                                         int offset_i, int offset_j)
    : MatrixView(mat, offset_i, offset_j, mat.nrRows() - offset_i, mat.nrCols() - offset_j) {
  assert(offset_i < mat.nrCols());
  assert(offset_j < mat.nrRows());
}

template <class Scalar, DeviceType device_t>
template <template <class, DeviceType> class Matrix>
MatrixView<Scalar, device_t>::MatrixView(Matrix<Scalar, device_t>& mat, int offset_i, int offset_j,
                                         int ni, int nj)
    : ptr_(mat.ptr(offset_i, offset_j)), ldm_(mat.leadingDimension()), size_(std::make_pair(ni, nj)) {
  assert(ni + offset_i <= mat.nrRows());
  assert(nj + offset_j <= mat.nrCols());
}
template <class Scalar, DeviceType device_t>
template <template <class, DeviceType> class Matrix>
MatrixView<Scalar, device_t>::MatrixView(const Matrix<std::remove_cv_t<Scalar>, device_t>& mat,
                                         int offset_i, int offset_j, int ni, int nj)
    : ptr_(mat.ptr(offset_i, offset_j)), ldm_(mat.leadingDimension()), size_(std::make_pair(ni, nj)) {
  assert(ni + offset_i <= mat.nrRows());
  assert(nj + offset_j <= mat.nrCols());
}

template <class Scalar, DeviceType device_t>
template <template <class, DeviceType> class Matrix, class Scalar2>
MatrixView<Scalar, device_t>& MatrixView<Scalar, device_t>::operator=(
    const Matrix<Scalar2, device_t>& rhs) {
  static_assert(device_t == CPU, "Copy implemented only on CPU");
  if (nrCols() != rhs.nrCols() || nrRows() != rhs.nrRows()) {
    throw(std::invalid_argument("Matrix size mismatch."));
  }

  for (int j = 0; j < nrCols(); ++j)
    for (int i = 0; i < nrRows(); ++i)
      (*this)(i, j) = rhs(i, j);
  return *this;
}

template <class Scalar, DeviceType device_t>
MatrixView<Scalar, device_t>& MatrixView<Scalar, device_t>::operator=(const MatrixView& rhs) {
  static_assert(device_t == CPU, "Copy implemented only on CPU");
  if (nrCols() != rhs.nrCols() || nrRows() != rhs.nrRows()) {
    throw(std::invalid_argument("Matrix size mismatch."));
  }

  for (int j = 0; j < nrCols(); ++j)
    for (int i = 0; i < nrRows(); ++i)
      (*this)(i, j) = rhs(i, j);
  return *this;
}

template <class Scalar, DeviceType device_t>
void MatrixView<Scalar, device_t>::print(std::ostream& out) const {
  out << "\tMatrix view:\n";
  out << "Size: \t" << size_.first << ", " << size_.second << "\n";

  for (int j = 0; j < nrRows(); ++j) {
    for (int i = 0; i < nrRows(); ++i)
      out << (*this)(i, j) << "\t";
    out << "\n";
  }
  out << "\n" << std::endl;
}

}  // namespace linalg
}  // namespace dca

#endif  // DCA_LINALG_MATRIX_VIEW_HPP
