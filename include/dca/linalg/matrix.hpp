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
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include "dca/linalg/device_type.hpp"

#include "comp_library/linalg/src/linalg_operations/copy_from_tem.h"
#include "comp_library/linalg/src/linalg_operations/copy_from_CPU_CPU.h"
#include "comp_library/linalg/src/linalg_operations/copy_from_CPU_GPU.h"
#include "comp_library/linalg/src/linalg_operations/copy_from_GPU_CPU.h"
#include "comp_library/linalg/src/linalg_operations/copy_from_GPU_GPU.h"
#include "comp_library/linalg/src/matrix_scalartype.h"
#include "comp_library/linalg/src/linalg_structures/cublas_thread_manager_tem.h"
#include "comp_library/linalg/src/linalg_structures/cublas_thread_manager_CPU.h"
#include "comp_library/linalg/src/linalg_structures/cublas_thread_manager_GPU.h"
#include "comp_library/linalg/src/linalg_operations/memory_management_CPU.h"
#include "comp_library/linalg/src/linalg_operations/memory_management_GPU.h"

namespace dca {
namespace linalg {
// dca::linalg::

using namespace ::LIN_ALG;

template <typename ScalarType, DeviceType device_name>
class Matrix {
public:
  using ThisType = Matrix<ScalarType, device_name>;
  using ValueType = typename MATRIX_SCALARTYPE<ScalarType, device_name>::new_scalartype;

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
  Matrix(std::string name, std::pair<int, int> size, std::pair<int, int> capacity, int thread_id,
         int stream_id);

  Matrix(const Matrix<ScalarType, device_name>& rhs);

  template <DeviceType rhs_device_name>
  Matrix(const Matrix<ScalarType, rhs_device_name>& rhs);

  ~Matrix();

  Matrix<ScalarType, device_name>& operator=(const Matrix<ScalarType, device_name>& rhs);

  template <DeviceType rhs_device_name>
  Matrix<ScalarType, device_name>& operator=(const Matrix<ScalarType, rhs_device_name>& rhs);

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

  // TODO: remove reference (needed to check for external changes)
  const int& get_thread_id() const {
    return thread_id_;
  }
  const int& get_stream_id() const {
    return stream_id_;
  }

  void setThreadAndStreamId(int thread_id, int stream_id) {
    thread_id_ = thread_id;
    stream_id_ = stream_id;
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

  // Copies the value of rhs inside *this.
  // rhs and *this must have the same size.
  template <DeviceType rhs_device_name>
  void copy_from(Matrix<ScalarType, rhs_device_name>& rhs);

  // Copies the value of rhs inside *this.
  // rhs and *this must have the same size.
  template <DeviceType rhs_device_name>
  void copy_from(Matrix<ScalarType, rhs_device_name>& rhs, copy_concurrency_type copy_t);

  // Prints the values of the matrix elements.
  void print() const;
  // Prints the properties of *this.
  void print_fingerprint() const;

  // TODO: move to matrix operations
  template <DeviceType rhs_device_name>
  ScalarType difference(Matrix<ScalarType, rhs_device_name>& rhs, double diff_threshold = 1e-3) const;

private:
  static std::pair<int, int> capacityMultipleOfBlockSize(std::pair<int, int> size);
  inline static size_t nrElements(std::pair<int, int> size) {
    return static_cast<size_t>(size.first) * static_cast<size_t>(size.second);
  }
  const static int block_size = 32;

  std::string name_;

  std::pair<int, int> size_;
  std::pair<int, int> capacity_;

  int thread_id_;
  int stream_id_;

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
    : Matrix(str, size, capacity, -1, -1) {}

template <typename ScalarType, DeviceType device_name>
Matrix<ScalarType, device_name>::Matrix(std::string str, std::pair<int, int> size,
                                        std::pair<int, int> capacity, int thread_id, int stream_id)
    : name_(str),
      size_(size),
      capacity_(capacityMultipleOfBlockSize(capacity)),
      thread_id_(thread_id),
      stream_id_(stream_id),
      data_(nullptr) {
  assert(size_.first >= 0 && size_.second >= 0);
  assert(capacity.first >= 0 && capacity.second >= 0);
  assert(capacity.first >= size_.first && capacity.second >= size_.second);
  assert(capacity_.first >= capacity.first && capacity_.second >= capacity.second);

  MEMORY_MANAGEMENT<device_name>::allocate(data_, capacity_);
  MEMORY_MANAGEMENT<device_name>::set_to_zero(data_, nrElements(capacity_));
}

template <typename ScalarType, DeviceType device_name>
Matrix<ScalarType, device_name>::Matrix(const Matrix<ScalarType, device_name>& rhs)
    : name_(rhs.name_),
      size_(rhs.size_),
      capacity_(rhs.capacity_),
      thread_id_(-1),
      stream_id_(-1),
      data_(nullptr) {
  MEMORY_MANAGEMENT<device_name>::allocate(data_, capacity_);

  COPY_FROM<device_name, device_name>::execute(rhs.data_, rhs.size_, rhs.capacity_, data_, size_,
                                               capacity_);
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType rhs_device_name>
Matrix<ScalarType, device_name>::Matrix(const Matrix<ScalarType, rhs_device_name>& rhs)
    : name_(rhs.name_),
      size_(rhs.size_),
      capacity_(rhs.capacity_),
      thread_id_(-1),
      stream_id_(-1),
      data_(nullptr) {
  MEMORY_MANAGEMENT<device_name>::allocate(data_, capacity_);

  COPY_FROM<rhs_device_name, device_name>::execute(rhs.data_, rhs.size_, rhs.capacity_, data_,
                                                   size_, capacity_);
}

template <typename ScalarType, DeviceType device_name>
Matrix<ScalarType, device_name>::~Matrix() {
  MEMORY_MANAGEMENT<device_name>::deallocate(data_);
}

template <typename ScalarType, DeviceType device_name>
void Matrix<ScalarType, device_name>::resize(std::pair<int, int> new_size) {
  assert(new_size.first >= 0 && new_size.second >= 0);
  if (new_size.first > capacity_.first || new_size.second > capacity_.second) {
    // CUBLAS_THREAD_MANAGER<device_name>::synchronize_streams(thread_id_, stream_id_);
    std::pair<int, int> new_capacity = capacityMultipleOfBlockSize(new_size);

    ValueType* new_data = NULL;
    MEMORY_MANAGEMENT<device_name>::allocate(new_data, new_capacity);
    COPY_FROM<device_name, device_name>::execute(data_, size_, capacity_, new_data, size_,
                                                 new_capacity);
    MEMORY_MANAGEMENT<device_name>::deallocate(data_);

    data_ = new_data;
    capacity_ = new_capacity;
    size_ = new_size;

    // CUBLAS_THREAD_MANAGER<device_name>::synchronize_streams(thread_id_, stream_id_);
  }
  else {
    size_ = new_size;
  }
}

template <typename ScalarType, DeviceType device_name>
Matrix<ScalarType, device_name>& Matrix<ScalarType, device_name>::operator=(
    const Matrix<ScalarType, device_name>& rhs) {
  name_ = rhs.name_;
  resizeNoCopy(rhs.size_);
  thread_id_ = -1;
  stream_id_ = -1;

  COPY_FROM<device_name, device_name>::execute(rhs.data_, rhs.size_, rhs.capacity_, data_, size_,
                                               capacity_);
  return *this;
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType rhs_device_name>
Matrix<ScalarType, device_name>& Matrix<ScalarType, device_name>::operator=(
    const Matrix<ScalarType, rhs_device_name>& rhs) {
  name_ = rhs.name_;
  resizeNoCopy(rhs.size_);
  thread_id_ = -1;
  stream_id_ = -1;

  COPY_FROM<rhs_device_name, device_name>::execute(rhs.data_, rhs.size_, rhs.capacity_, data_,
                                                   size_, capacity_);
  return *this;
}

template <typename ScalarType, DeviceType device_name>
void Matrix<ScalarType, device_name>::resizeNoCopy(std::pair<int, int> new_size) {
  if (new_size.first > capacity_.first || new_size.second > capacity_.second) {
    // CUBLAS_THREAD_MANAGER<device_name>::synchronize_streams(thread_id_, stream_id_);

    size_ = new_size;
    capacity_ = capacityMultipleOfBlockSize(new_size);

    MEMORY_MANAGEMENT<device_name>::deallocate(data_);
    MEMORY_MANAGEMENT<device_name>::allocate(data_, capacity_);

    // CUBLAS_THREAD_MANAGER<device_name>::synchronize_streams(thread_id_, stream_id_);
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
  std::swap(thread_id_, rhs.thread_id_);
  std::swap(stream_id_, rhs.stream_id_);
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType rhs_device_name>
void Matrix<ScalarType, device_name>::copy_from(Matrix<ScalarType, rhs_device_name>& rhs) {
  resize(rhs.size_);

  COPY_FROM<rhs_device_name, device_name>::execute(rhs.data_, rhs.size_, rhs.capacity_, data_,
                                                   size_, capacity_);
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType rhs_device_name>
void Matrix<ScalarType, device_name>::copy_from(Matrix<ScalarType, rhs_device_name>& rhs,
                                                copy_concurrency_type copy_t) {
  const static DeviceType device_t =
      LIN_ALG::CUBLAS_DEVICE_NAME<rhs_device_name, device_name>::device_t;

  assert(thread_id_ > -1 and stream_id_ > -1);

  resize(rhs.size_);

  switch (copy_t) {
    case SYNCHRONOUS:
      COPY_FROM<rhs_device_name, device_name>::execute(rhs.data_, rhs.size_, rhs.capacity_, data_,
                                                       size_, capacity_);
      break;

    case ASYNCHRONOUS:
      CUBLAS_THREAD_MANAGER<device_t>::synchronize_streams(thread_id_, stream_id_);

      COPY_FROM<rhs_device_name, device_name>::execute(rhs.data_, rhs.size_, rhs.capacity_, data_,
                                                       size_, capacity_, thread_id_, stream_id_);

      CUBLAS_THREAD_MANAGER<device_t>::synchronize_streams(thread_id_, stream_id_);
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType rhs_device_name>
ScalarType Matrix<ScalarType, device_name>::difference(Matrix<ScalarType, rhs_device_name>& rhs,
                                                       double diff_threshold) const {
  if (size_ != rhs.size_) {
    throw std::logic_error("different matrix size");
  }

  Matrix<ScalarType, CPU> cp_this(*this);
  Matrix<ScalarType, CPU> cp_rhs(rhs);

  auto max_diff = std::abs(ScalarType(0));

  for (int j = 0; j < size_.second; ++j) {
    for (int i = 0; i < size_.first; ++i) {
      max_diff = std::max(max_diff, std::fabs(cp_this(i, j) - cp_rhs(i, j)));
    }
  }

  if (std::fabs(max_diff) > diff_threshold) {
#ifndef DNDEBUG
    std::stringstream s;
    for (int i = 0; i < size_.first; ++i) {
      for (int j = 0; j < size_.second; ++j) {
        if (std::fabs(cp_this(i, j) - cp_rhs(i, j)) <= diff_threshold)
          s << 0. << "\t";
        else
          s << cp_this(i, j) - cp_rhs(i, j) << "\t";
      }
      s << "\n";
    }
    std::cout << s.str() << std::endl;
#endif  // DNDEBUG

    throw std::logic_error(__FUNCTION__);
  }

  return max_diff;
}

template <typename ScalarType, DeviceType device_name>
void Matrix<ScalarType, device_name>::print() const {
  MEMORY_MANAGEMENT<device_name>::print(data_, size_, capacity_);
}

template <typename ScalarType, DeviceType device_name>
void Matrix<ScalarType, device_name>::print_fingerprint() const {
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

  size.first = (size.first + block_size - 1) / block_size * block_size;
  size.second = (size.second + block_size - 1) / block_size * block_size;

  return size;
}

}  // linalg
}  // dca

namespace LIN_ALG {
template <typename ScalarType, dca::linalg::DeviceType device_name>
using matrix = typename dca::linalg::Matrix<ScalarType, device_name>;
}

#endif  // DCA_LINALG_MATRIX_HPP
