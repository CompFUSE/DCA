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
// This file provides the Matrix object for different device types.

#ifndef DCA_LINALG_MATRIX_HPP
#define DCA_LINALG_MATRIX_HPP

#include <cassert>
#include <cmath>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include "dca/linalg/device_type.hpp"
#include "dca/linalg/util/copy.hpp"
#include "dca/linalg/util/memory.hpp"

namespace dca {
namespace linalg {
// dca::linalg::

template <typename ScalarType, DeviceType device_name>
class Matrix {
public:
  using ThisType = Matrix<ScalarType, device_name>;
  using ValueType = ScalarType;

  Matrix();
  Matrix(std::string name);

  Matrix(int size);
  Matrix(std::string name, int size);

  // Preconditions: capacity >= size.
  Matrix(int size, int capacity);
  Matrix(std::string name, int size, int capacity);

  Matrix(std::pair<int, int> size);
  Matrix(std::string name, std::pair<int, int> size);

  // Preconditions: capacity.first >= size.first, capacity.second >= size.second.
  Matrix(std::pair<int, int> size, std::pair<int, int> capacity);
  Matrix(std::string name, std::pair<int, int> size, std::pair<int, int> capacity);

  Matrix(const Matrix<ScalarType, device_name>& rhs);

  template <DeviceType rhs_device_name>
  Matrix(const Matrix<ScalarType, rhs_device_name>& rhs);

  ~Matrix();

  Matrix<ScalarType, device_name>& operator=(const Matrix<ScalarType, device_name>& rhs);

  template <DeviceType rhs_device_name>
  Matrix<ScalarType, device_name>& operator=(const Matrix<ScalarType, rhs_device_name>& rhs);

  // Returns true if this is equal to other, false otherwise.
  // Two matrices are equal, if they have the same size and contain the same elements. Name and
  // capacity are ignored.
  // Special case: two matrices without elements are equal.
  template <DeviceType dn = device_name>
  bool operator==(
      std::enable_if_t<device_name == CPU and dn == CPU, const Matrix<ScalarType, dn>&> other) const;
  template <DeviceType dn = device_name>
  // Returns true if this is not equal to other, false otherwise.
  // See description of operator== for the definition of equality.
  bool operator!=(
      std::enable_if_t<device_name == CPU and dn == CPU, const Matrix<ScalarType, dn>&> other) const;

  // Returns the (i,j)-th element of the matrix.
  // Preconditions: 0 <= i < size().first, 0 <= j < size().second.
  // This method is available only if device_name == CPU.
  template <DeviceType dn = device_name>
  std::enable_if_t<device_name == CPU && dn == CPU, ScalarType&> operator()(int i, int j) {
    assert(i >= 0 && i < size_.first);
    assert(j >= 0 && j < size_.second);
    return data_[i + j * leadingDimension()];
  }
  template <DeviceType dn = device_name>
  std::enable_if_t<device_name == CPU && dn == CPU, const ScalarType&> operator()(int i, int j) const {
    assert(i >= 0 && i < size_.first);
    assert(j >= 0 && j < size_.second);
    return data_[i + j * leadingDimension()];
  }

  const std::string& get_name() const {
    return name_;
  }
  void set_name(const std::string& new_name) {
    name_ = new_name;
  }

  // Returns the pointer to the (0,0)-th element.
  ValueType* ptr() {
    return data_;
  }
  const ValueType* ptr() const {
    return data_;
  }

  // Returns the pointer to the (i,j)-th element.
  // Preconditions: 0 <= i < size().first, 0 <= j < size().second.
  ValueType* ptr(int i, int j) {
    assert(i >= 0 && i < size_.first);
    assert(j >= 0 && j < size_.second);
    return data_ + i + j * leadingDimension();
  }
  const ValueType* ptr(int i, int j) const {
    assert(i >= 0 && i < size_.first);
    assert(j >= 0 && j < size_.second);
    return data_ + i + j * leadingDimension();
  }

  bool is_square() const {
    return (size_.first == size_.second);
  }

  // TODO: remove reference (needed to check for external changes)
  const std::pair<int, int>& size() const {
    return size_;
  }
  const std::pair<int, int>& capacity() const {
    return capacity_;
  }
  int nrRows() const {
    return size_.first;
  }
  int nrCols() const {
    return size_.second;
  }
  int leadingDimension() const {
    return capacity_.first;
  }

  // Resizes *this to a (new_size * new_size) matrix.
  // Elements added may have any value.
  // Remark: The capacity of the matrix and element pointers do not change
  // if new_size <= capacity().first and new_size <= capacity().second.
  void resize(int new_size) {
    resize(std::make_pair(new_size, new_size));
  }
  // Resizes *this to a (new_size.first * new_size.second) matrix.
  // Elements added may have any value.
  // Remark: The capacity of the matrix and element pointers do not change
  // if new_size.first <= capacity().first and new_size.second <= capacity().second.
  void resize(std::pair<int, int> new_size);

  // Resizes *this to a (new_size * new_size) matrix.
  // The previous elements are not copied, therefore all the elements
  // may have any value after the call to this method.
  // Remark: The capacity of the matrix and element pointers do not change
  // if new_size <= capacity().first and new_size <= capacity().second.
  void resizeNoCopy(int new_size) {
    resizeNoCopy(std::make_pair(new_size, new_size));
  }
  // Resizes *this to a (new_size.first * new_size.second) matrix.
  // The previous elements are not copied, therefore all the elements
  // may have any value after the call to this method.
  // Remark: The capacity of the matrix and element pointers do not change
  // if new_size.first <= capacity().first and new_size.second <= capacity().second.
  void resizeNoCopy(std::pair<int, int> new_size);

  void swap(Matrix<ScalarType, device_name>& rhs);

  // Asynchronous assignement (copy with stream = getStream(thread_id, stream_id))
  // + synchronization of stream
  // Preconditions: 0 <= thread_id < DCA_MAX_THREADS,
  //                0 <= stream_id < DCA_STREAMS_PER_THREADS.
  template <DeviceType rhs_device_name>
  void set(const Matrix<ScalarType, rhs_device_name>& rhs, int thread_id, int stream_id);

  // Prints the values of the matrix elements.
  template <DeviceType dn = device_name>
  std::enable_if_t<device_name == CPU && dn == CPU, void> print() const;
  template <DeviceType dn = device_name>
  std::enable_if_t<device_name != CPU && dn == device_name, void> print() const;
  // Prints the properties of *this.
  void printFingerprint() const;

private:
  static std::pair<int, int> capacityMultipleOfBlockSize(std::pair<int, int> size);
  inline static size_t nrElements(std::pair<int, int> size) {
    return static_cast<size_t>(size.first) * static_cast<size_t>(size.second);
  }
  const static int block_size_ = 32;

  std::string name_;

  std::pair<int, int> size_;
  std::pair<int, int> capacity_;

  ValueType* data_;

  template <class ScalarType2, DeviceType device_name2>
  friend class dca::linalg::Matrix;
};

template <typename ScalarType, DeviceType device_name>
Matrix<ScalarType, device_name>::Matrix() : Matrix(0) {}

template <typename ScalarType, DeviceType device_name>
Matrix<ScalarType, device_name>::Matrix(std::string str) : Matrix(str, 0) {}

template <typename ScalarType, DeviceType device_name>
Matrix<ScalarType, device_name>::Matrix(int size) : Matrix(std::make_pair(size, size)) {}

template <typename ScalarType, DeviceType device_name>
Matrix<ScalarType, device_name>::Matrix(std::string str, int size)
    : Matrix(str, std::make_pair(size, size)) {}

template <typename ScalarType, DeviceType device_name>
Matrix<ScalarType, device_name>::Matrix(int size, int capacity)
    : Matrix(std::make_pair(size, size), std::make_pair(capacity, capacity)) {}

template <typename ScalarType, DeviceType device_name>
Matrix<ScalarType, device_name>::Matrix(std::string str, int size, int capacity)
    : Matrix(str, std::make_pair(size, size), std::make_pair(capacity, capacity)) {}

template <typename ScalarType, DeviceType device_name>
Matrix<ScalarType, device_name>::Matrix(std::pair<int, int> size) : Matrix(size, size) {}

template <typename ScalarType, DeviceType device_name>
Matrix<ScalarType, device_name>::Matrix(std::string str, std::pair<int, int> size)
    : Matrix(str, size, size) {}

template <typename ScalarType, DeviceType device_name>
Matrix<ScalarType, device_name>::Matrix(std::pair<int, int> size, std::pair<int, int> capacity)
    : Matrix("unnamed matrix", size, capacity) {}

template <typename ScalarType, DeviceType device_name>
Matrix<ScalarType, device_name>::Matrix(std::string str, std::pair<int, int> size,
                                        std::pair<int, int> capacity)
    : name_(str), size_(size), capacity_(capacityMultipleOfBlockSize(capacity)), data_(nullptr) {
  assert(size_.first >= 0 && size_.second >= 0);
  assert(capacity.first >= 0 && capacity.second >= 0);
  assert(capacity.first >= size_.first && capacity.second >= size_.second);
  assert(capacity_.first >= capacity.first && capacity_.second >= capacity.second);

  util::Memory<device_name>::allocate(data_, nrElements(capacity_));
  util::Memory<device_name>::setToZero(data_, nrElements(capacity_));
}

template <typename ScalarType, DeviceType device_name>
Matrix<ScalarType, device_name>::Matrix(const Matrix<ScalarType, device_name>& rhs)
    : name_(rhs.name_), size_(rhs.size_), capacity_(rhs.capacity_), data_(nullptr) {
  util::Memory<device_name>::allocate(data_, nrElements(capacity_));
  util::memoryCopy(data_, leadingDimension(), rhs.data_, rhs.leadingDimension(), size_);
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType rhs_device_name>
Matrix<ScalarType, device_name>::Matrix(const Matrix<ScalarType, rhs_device_name>& rhs)
    : name_(rhs.name_), size_(rhs.size_), capacity_(rhs.capacity_), data_(nullptr) {
  util::Memory<device_name>::allocate(data_, nrElements(capacity_));
  util::memoryCopy(data_, leadingDimension(), rhs.data_, rhs.leadingDimension(), size_);
}

template <typename ScalarType, DeviceType device_name>
Matrix<ScalarType, device_name>::~Matrix() {
  util::Memory<device_name>::deallocate(data_);
}

template <typename ScalarType, DeviceType device_name>
void Matrix<ScalarType, device_name>::resize(std::pair<int, int> new_size) {
  assert(new_size.first >= 0 && new_size.second >= 0);
  if (new_size.first > capacity_.first || new_size.second > capacity_.second) {
    std::pair<int, int> new_capacity = capacityMultipleOfBlockSize(new_size);

    ValueType* new_data = NULL;
    util::Memory<device_name>::allocate(new_data, nrElements(new_capacity));
    util::memoryCopy(new_data, new_capacity.first, data_, leadingDimension(), size_);
    util::Memory<device_name>::deallocate(data_);

    data_ = new_data;
    capacity_ = new_capacity;
    size_ = new_size;
  }
  else {
    size_ = new_size;
  }
}

template <typename ScalarType, DeviceType device_name>
Matrix<ScalarType, device_name>& Matrix<ScalarType, device_name>::operator=(
    const Matrix<ScalarType, device_name>& rhs) {
  resizeNoCopy(rhs.size_);
  util::memoryCopy(data_, leadingDimension(), rhs.data_, rhs.leadingDimension(), size_);
  return *this;
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType rhs_device_name>
Matrix<ScalarType, device_name>& Matrix<ScalarType, device_name>::operator=(
    const Matrix<ScalarType, rhs_device_name>& rhs) {
  resizeNoCopy(rhs.size_);
  util::memoryCopy(data_, leadingDimension(), rhs.data_, rhs.leadingDimension(), size_);
  return *this;
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType dn>
bool Matrix<ScalarType, device_name>::operator==(
    std::enable_if_t<device_name == CPU and dn == CPU, const Matrix<ScalarType, dn>&> other) const {
  if (size() != other.size())
    return nrRows() * nrCols() == 0 and other.nrRows() * other.nrCols() == 0;

  for (int j = 0; j < nrCols(); ++j)
    for (int i = 0; i < nrRows(); ++i)
      if ((*this)(i, j) != other(i, j))
        return false;

  return true;
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType dn>
bool Matrix<ScalarType, device_name>::operator!=(
    std::enable_if_t<device_name == CPU and dn == CPU, const Matrix<ScalarType, dn>&> other) const {
  return not(*this == other);
}

template <typename ScalarType, DeviceType device_name>
void Matrix<ScalarType, device_name>::resizeNoCopy(std::pair<int, int> new_size) {
  if (new_size.first > capacity_.first || new_size.second > capacity_.second) {
    size_ = new_size;
    capacity_ = capacityMultipleOfBlockSize(new_size);

    util::Memory<device_name>::deallocate(data_);
    util::Memory<device_name>::allocate(data_, nrElements(capacity_));
  }
  else {
    size_ = new_size;
  }
}

template <typename ScalarType, DeviceType device_name>
void Matrix<ScalarType, device_name>::swap(Matrix<ScalarType, device_name>& rhs) {
  std::swap(name_, rhs.name_);
  std::swap(size_, rhs.size_);
  std::swap(capacity_, rhs.capacity_);
  std::swap(data_, rhs.data_);
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType rhs_device_name>
void Matrix<ScalarType, device_name>::set(const Matrix<ScalarType, rhs_device_name>& rhs,
                                          int thread_id, int stream_id) {
  resize(rhs.size_);
  util::memoryCopy(data_, leadingDimension(), rhs.data_, rhs.leadingDimension(), size_, thread_id,
                   stream_id);
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType dn>
std::enable_if_t<device_name == CPU && dn == CPU, void> Matrix<ScalarType, device_name>::print() const {
  printFingerprint();

  std::stringstream ss;
  ss.precision(6);
  ss << std::scientific;

  ss << "\n";
  for (int i = 0; i < nrRows(); ++i) {
    for (int j = 0; j < nrCols(); ++j)
      ss << "\t" << operator()(i, j);
    ss << "\n";
  }

  std::cout << ss.str() << std::endl;
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType dn>
std::enable_if_t<device_name != CPU && dn == device_name, void> Matrix<ScalarType,
                                                                       device_name>::print() const {
  Matrix<ScalarType, CPU> copy(*this);
  copy.print();
}

template <typename ScalarType, DeviceType device_name>
void Matrix<ScalarType, device_name>::printFingerprint() const {
  std::stringstream ss;

  ss << "\n";
  ss << "    name: " << name_ << "\n";
  ss << "    size: " << size_.first << ", " << size_.second << "\n";
  ss << "    capacity: " << capacity_.first << ", " << capacity_.second << "\n";
  ss << "    memory-size: " << nrElements(capacity_) * sizeof(ScalarType) * 1.e-6 << "(Mbytes)\n";

  std::cout << ss.str() << std::endl;
}

template <typename ScalarType, DeviceType device_name>
std::pair<int, int> Matrix<ScalarType, device_name>::capacityMultipleOfBlockSize(
    std::pair<int, int> size) {
  assert(size.first >= 0);
  assert(size.second >= 0);

  size.first = (size.first + block_size_ - 1) / block_size_ * block_size_;
  size.second = (size.second + block_size_ - 1) / block_size_ * block_size_;

  return size;
}

}  // linalg
}  // dca

#endif  // DCA_LINALG_MATRIX_HPP
