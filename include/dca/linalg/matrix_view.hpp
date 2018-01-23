// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi(gbalduzz@itp.phys.ethz.ch)
//
// This file provides a lightweight proxy to access blocks of a matrix. The underlying matrix must
// not be destroyed while this object is in use. A const and a non const version are provided.

#ifndef DCA_LINALG_MATRIX_VIEW_HPP
#define DCA_LINALG_MATRIX_VIEW_HPP

#include <cassert>
#include <iostream>
#include <stdexcept>

#include "dca/linalg/device_type.hpp"
#include "dca/linalg/matrix.hpp"

namespace dca {
namespace linalg {
// dca::linalg::

template <typename ScalarType, DeviceType device_name = linalg::CPU>
class MatrixView {
public:
  using Pair = std::pair<int, int>;
  MatrixView(ScalarType* data, Pair size, int ld);
  MatrixView(ScalarType* data, int size, int ld);
  MatrixView(ScalarType* data, int size);
  MatrixView(Matrix<ScalarType, device_name>& mat);
  MatrixView(Matrix<ScalarType, device_name>& mat, int offset_i, int offset_j);
  MatrixView(Matrix<ScalarType, device_name>& mat, int offset_i, int offset_j, int ni, int nj);

  template <class MatrixB>
  MatrixView& operator=(const MatrixB& rhs);

  void print(std::ostream& out = std::cout) const;

  inline int leadingDimension() const {
    return ldm_;
  }
  inline int ld() const {
    return leadingDimension();
  }

  std::pair<int, int> size() const {
    return size_;
  }
  ScalarType* ptr() {
    return ptr_;
  }
  const ScalarType* ptr() const {
    return ptr_;
  }
  ScalarType* ptr(int i, int j) {
    return ptr_ + leadingDimension() * j + i;
  }
  const ScalarType* ptr(int i, int j) const {
    return ptr_ + leadingDimension() * j + i;
  }
  int nrRows() const {
    return size_.first;
  }
  int nrCols() const {
    return size_.second;
  }
  bool is_square() const {
    return size_.first == size_.second;
  }
  ScalarType& operator()(int i, int j) {
    return ptr_[i + j * ldm_];
  }
  const ScalarType& operator()(int i, int j) const {
    return ptr_[i + j * ldm_];
  }

private:
  ScalarType* const ptr_;
  const int ldm_;
  const std::pair<int, int> size_;
};

template <typename ScalarType, DeviceType device_name = linalg::CPU>
class ConstMatrixView {
public:
  using Pair = std::pair<int, int>;
  ConstMatrixView(const ScalarType* data, Pair size, int ld);
  ConstMatrixView(const ScalarType* data, int size, int ld);
  ConstMatrixView(const ScalarType* data, int size);
  ConstMatrixView(const Matrix<ScalarType, device_name>& mat);
  ConstMatrixView(const Matrix<ScalarType, device_name>& mat, int offset_i, int offset_j);
  ConstMatrixView(const Matrix<ScalarType, device_name>& mat, int offset_i, int offset_j, int ni,
                  int nj);

  void print(std::ostream& out = std::cout) const;

  inline int leadingDimension() const {
    return ldm_;
  }
  inline int ld() const {
    return leadingDimension();
  }

  std::pair<int, int> size() const {
    return size_;
  }
  ScalarType* ptr() {
    return ptr_;
  }
  const ScalarType* ptr() const {
    return ptr_;
  }
  const ScalarType* ptr(int i, int j) const {
    return ptr_ + leadingDimension() * j + i;
  }
  int nrRows() const {
    return size_.first;
  }
  int nrCols() const {
    return size_.second;
  }
  bool is_square() const {
    return size_.first == size_.second;
  }

  const ScalarType& operator()(int i, int j) const {
    return ptr_[i + j * ldm_];
  }

private:
  const ScalarType* const ptr_;
  const int ldm_;
  const std::pair<int, int> size_;
};

// Implementation of non-const version.
template <typename ScalarType, DeviceType device_t>
MatrixView<ScalarType, device_t>::MatrixView(ScalarType* const data, const Pair size, const int ld)
    : ptr_(data), ldm_(ld), size_(size) {}

template <typename ScalarType, DeviceType device_t>
MatrixView<ScalarType, device_t>::MatrixView(ScalarType* const data, const int size, const int ld)
    : ptr_(data), ldm_(ld), size_(std::make_pair(size, size)) {}

template <typename ScalarType, DeviceType device_t>
MatrixView<ScalarType, device_t>::MatrixView(ScalarType* const data, const int size)
    : ptr_(data), ldm_(size), size_(std::make_pair(size, size)) {}

template <typename ScalarType, DeviceType device_t>
MatrixView<ScalarType, device_t>::MatrixView(Matrix<ScalarType, device_t>& mat)
    : ptr_(mat.ptr()), ldm_(mat.leadingDimension()), size_(mat.size()) {}

template <typename ScalarType, DeviceType device_t>
MatrixView<ScalarType, device_t>::MatrixView(Matrix<ScalarType, device_t>& mat, int offset_i,
                                             int offset_j)
    : ptr_(mat.ptr(offset_i, offset_j)), ldm_(mat.leadingDimension()), size_(mat.size()) {}

template <typename ScalarType, DeviceType device_t>
MatrixView<ScalarType, device_t>::MatrixView(Matrix<ScalarType, device_t>& mat, int offset_i,
                                             int offset_j, int ni, int nj)
    : ptr_(mat.ptr(offset_i, offset_j)),
      ldm_(mat.leadingDimension()),
      size_(std::make_pair(ni, nj)) {}

template <typename ScalarType, DeviceType device_t>
template <class MatrixB>
MatrixView<ScalarType, device_t>& MatrixView<ScalarType, device_t>::operator=(const MatrixB& rhs) {
  assert(nrCols() == rhs.nrCols() and nrRows() == rhs.nrRows());
  for (int j = 0; j < nrCols(); ++j)
    for (int i = 0; i < nrRows(); ++i)
      (*this)(i, j) = rhs(i, j);
  return *this;
}

template <typename ScalarType, DeviceType device_t>
void MatrixView<ScalarType, device_t>::print(std::ostream& out) const {
  out << "\tMatrix view:\n";
  out << "Size: \t" << size_.first << ", " << size_.second << "\n";

  for (int j = 0; j < nrRows(); ++j) {
    for (int i = 0; i < nrRows(); ++i)
      out << (*this)(i, j) << "\t";
    out << "\n";
  }
  out << "\n" << std::endl;
}

// Implementation of const version.

template <typename ScalarType, DeviceType device_t>
ConstMatrixView<ScalarType, device_t>::ConstMatrixView(const ScalarType* const data,
                                                       const Pair size, const int ld)
    : ptr_(data), ldm_(ld), size_(size) {}

template <typename ScalarType, DeviceType device_t>
ConstMatrixView<ScalarType, device_t>::ConstMatrixView(const ScalarType* const data, const int size,
                                                       const int ld)
    : ptr_(data), ldm_(ld), size_(std::make_pair(size, size)) {}

template <typename ScalarType, DeviceType device_t>
ConstMatrixView<ScalarType, device_t>::ConstMatrixView(const ScalarType* const data, const int size)
    : ptr_(data), ldm_(size), size_(std::make_pair(size, size)) {}

template <typename ScalarType, DeviceType device_t>

ConstMatrixView<ScalarType, device_t>::ConstMatrixView(const Matrix<ScalarType, device_t>& mat)
    : ptr_(mat.ptr()), ldm_(mat.leadingDimension()), size_(mat.size()) {}

template <typename ScalarType, DeviceType device_t>
ConstMatrixView<ScalarType, device_t>::ConstMatrixView(const Matrix<ScalarType, device_t>& mat,
                                                       int offset_i, int offset_j)
    : ptr_(mat.ptr(offset_i, offset_j)), ldm_(mat.leadingDimension()), size_(mat.size()) {}

template <typename ScalarType, DeviceType device_t>
ConstMatrixView<ScalarType, device_t>::ConstMatrixView(const Matrix<ScalarType, device_t>& mat,
                                                       int offset_i, int offset_j, int ni, int nj)
    : ptr_(mat.ptr(offset_i, offset_j)),
      ldm_(mat.leadingDimension()),
      size_(std::make_pair(ni, nj)) {}

template <typename ScalarType, DeviceType device_t>
void ConstMatrixView<ScalarType, device_t>::print(std::ostream& out) const {
  out << "\tMatrix view:\n";
  out << "Size: \t" << size_.first << ", " << size_.second << "\n";

  for (int j = 0; j < nrRows(); ++j) {
    for (int i = 0; i < nrRows(); ++i)
      out << (*this)(i, j) << "\t";
    out << "\n";
  }
  out << "\n" << std::endl;
}

}  // linalg
}  // dca

#endif  // DCA_LINALG_MATRIX_VIEW_HPP
