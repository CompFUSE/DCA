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
// This file provides the Vector object for different device types.

#ifndef DCA_LINALG_VECTOR_HPP
#define DCA_LINALG_VECTOR_HPP

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

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
class Vector {
public:
  using ThisType = Vector<ScalarType, device_name>;
  using ValueType = typename MATRIX_SCALARTYPE<ScalarType, device_name>::new_scalartype;

  Vector();
  Vector(std::string name);

  Vector(size_t size);
  Vector(std::string name, size_t size);

  // Preconditions: capacity >= size.
  Vector(size_t size, size_t capacity);
  Vector(std::string name, size_t size, size_t capacity);
  Vector(std::string name, size_t size, size_t capacity, int thread_id, int stream_id);

  Vector(const Vector<ScalarType, device_name>& rhs);

  template <DeviceType rhs_device_name>
  Vector(const Vector<ScalarType, rhs_device_name>& rhs);

  ~Vector();

  Vector<ScalarType, device_name>& operator=(const Vector<ScalarType, device_name>& rhs);

  template <DeviceType rhs_device_name>
  Vector<ScalarType, device_name>& operator=(const Vector<ScalarType, rhs_device_name>& rhs);

  // Returns the i-th element of the vector.
  // Preconditions: 0 <= i < size().first.
  // This method is available only if device_name == CPU.
  template <DeviceType dn = device_name>
  std::enable_if_t<device_name == CPU && dn == CPU, ScalarType&> operator[](size_t i) {
    assert(i < size_);
    return data_[i];
  }
  template <DeviceType dn = device_name>
  std::enable_if_t<device_name == CPU && dn == CPU, const ScalarType&> operator[](size_t i) const {
    assert(i < size_);
    return data_[i];
  }

  const std::string& get_name() const {
    return name_;
  }

  int get_thread_id() const {
    return thread_id_;
  }
  int get_stream_id() const {
    return stream_id_;
  }

  void setThreadAndStreamId(int thread_id, int stream_id) {
    thread_id_ = thread_id;
    stream_id_ = stream_id;
  }

  // TODO: set functions will be modified after the copy routines are updated.
  void set(std::vector<ScalarType>& rhs);
  void set(std::vector<ScalarType>& rhs, copy_concurrency_type copy_t);

  template <DeviceType rhs_device_name>
  void set(Vector<ScalarType, rhs_device_name>& rhs);

  template <DeviceType rhs_device_name>
  void set(Vector<ScalarType, rhs_device_name>& rhs, copy_concurrency_type copy_t);

  // Returns the pointer to the 0-th element of the vector.
  ValueType* ptr() {
    return data_;
  }
  const ValueType* ptr() const {
    return data_;
  }

  // Returns the pointer to the i-th element of the vector.
  // Preconditions: 0 <= i < size().first.
  ValueType* ptr(size_t i) {
    assert(i < size_);
    return data_ + i;
  }
  const ValueType* ptr(size_t i) const {
    assert(i < size_);
    return data_ + i;
  }

  size_t size() const {
    return size_;
  }
  size_t capacity() const {
    return capacity_;
  }

  // Resizes *this to a new_size vector.
  // Elements added may have any value.
  // Remark: The capacity of the vector and element pointers do not change
  // if new_size <= capacity().
  void resize(size_t new_size);
  // Resizes *this to a new_size vector.
  // The previous elements are not copied, therefore all the elements
  // may have any value after the call to this method.
  // Remark: The capacity of the vector and element pointers do not change
  // if new_size <= capacity().
  void resizeNoCopy(size_t new_size);

  // TODO: move to vector operations.
  template <DeviceType rhs_device_name>
  ScalarType difference(Vector<ScalarType, rhs_device_name>& rhs_Vector) const;

  // TODO: move to vector operations.
  ScalarType difference(std::vector<ScalarType>& rhs_Vector) const;

  // Prints the values of the vector elements.
  void print();

private:
  std::string name_;

  size_t size_;
  size_t capacity_;

  int thread_id_;
  int stream_id_;

  ValueType* data_;

  template <typename ScalarType2, DeviceType device_name2>
  friend class dca::linalg::Vector;
};

template <typename ScalarType, DeviceType device_name>
Vector<ScalarType, device_name>::Vector() : Vector("unnamed Vector", 0, 64, -1, -1) {}

template <typename ScalarType, DeviceType device_name>
Vector<ScalarType, device_name>::Vector(std::string name) : Vector(name, 0, 64, -1, -1) {}

template <typename ScalarType, DeviceType device_name>
Vector<ScalarType, device_name>::Vector(size_t size)
    : Vector("unnamed Vector", size, size, -1, -1) {}

template <typename ScalarType, DeviceType device_name>
Vector<ScalarType, device_name>::Vector(std::string name, size_t size)
    : Vector(name, size, size, -1, -1) {}

template <typename ScalarType, DeviceType device_name>
Vector<ScalarType, device_name>::Vector(size_t size, size_t capacity)
    : Vector("unnamed Vector", size, capacity, -1, -1) {}

template <typename ScalarType, DeviceType device_name>
Vector<ScalarType, device_name>::Vector(std::string name, size_t size, size_t capacity)
    : Vector(name, size, capacity, -1, -1) {}

template <typename ScalarType, DeviceType device_name>
Vector<ScalarType, device_name>::Vector(std::string name, size_t size, size_t capacity,
                                        int thread_id, int stream_id)
    : name_(name),
      size_(size),
      capacity_(capacity),
      thread_id_(thread_id),
      stream_id_(stream_id),
      data_(nullptr) {
  assert(capacity_ >= size_);
  MEMORY_MANAGEMENT<device_name>::allocate(data_, capacity_);
  MEMORY_MANAGEMENT<device_name>::set_to_zero(data_, capacity_);
}

template <typename ScalarType, DeviceType device_name>
Vector<ScalarType, device_name>::Vector(const Vector<ScalarType, device_name>& rhs)
    : name_(rhs.name_),
      size_(rhs.size_),
      capacity_(rhs.capacity_),
      thread_id_(-1),
      stream_id_(-1),
      data_(nullptr) {
  assert(capacity_ >= size_);
  MEMORY_MANAGEMENT<device_name>::allocate(data_, capacity_);
  COPY_FROM<device_name, device_name>::execute(rhs.data_, data_, size_);
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType rhs_device_name>
Vector<ScalarType, device_name>::Vector(const Vector<ScalarType, rhs_device_name>& rhs)
    : name_(rhs.name_),
      size_(rhs.size_),
      capacity_(rhs.capacity_),
      thread_id_(-1),
      stream_id_(-1),
      data_(nullptr) {
  assert(capacity_ >= size_);
  MEMORY_MANAGEMENT<device_name>::allocate(data_, capacity_);
  COPY_FROM<rhs_device_name, device_name>::execute(rhs.data_, data_, size_);
}

template <typename ScalarType, DeviceType device_name>
Vector<ScalarType, device_name>::~Vector() {
  MEMORY_MANAGEMENT<device_name>::deallocate(data_);
}

template <typename ScalarType, DeviceType device_name>
Vector<ScalarType, device_name>& Vector<ScalarType, device_name>::operator=(
    const Vector<ScalarType, device_name>& rhs) {
  name_ = rhs.name_;
  resizeNoCopy(rhs.size_);
  thread_id_ = -1;
  stream_id_ = -1;

  COPY_FROM<device_name, device_name>::execute(rhs.data_, data_, size_);
  return *this;
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType rhs_device_name>
Vector<ScalarType, device_name>& Vector<ScalarType, device_name>::operator=(
    const Vector<ScalarType, rhs_device_name>& rhs) {
  name_ = rhs.name_;
  resizeNoCopy(rhs.size_);
  thread_id_ = -1;
  stream_id_ = -1;

  COPY_FROM<rhs_device_name, device_name>::execute(rhs.data_, data_, size_);
  return *this;
}

template <typename ScalarType, DeviceType device_name>
void Vector<ScalarType, device_name>::set(std::vector<ScalarType>& rhs) {
  resizeNoCopy(rhs.size());
  COPY_FROM<dca::linalg::CPU, device_name>::execute(&rhs[0], data_, size_);
}

template <typename ScalarType, DeviceType device_name>
void Vector<ScalarType, device_name>::set(std::vector<ScalarType>& rhs, copy_concurrency_type copy_t) {
  const static DeviceType device_t =
      LIN_ALG::CUBLAS_DEVICE_NAME<dca::linalg::CPU, device_name>::device_t;

  resizeNoCopy(rhs.size());
  assert(thread_id_ > -1 and stream_id_ > -1);

  switch (copy_t) {
    case SYNCHRONOUS:
      COPY_FROM<dca::linalg::CPU, device_name>::execute(&rhs[0], data_, size_);
      break;

    case ASYNCHRONOUS:

      CUBLAS_THREAD_MANAGER<device_t>::synchronize_streams(thread_id_, stream_id_);

      COPY_FROM<dca::linalg::CPU, device_name>::execute(&rhs[0], data_, size_, thread_id_,
                                                        stream_id_);

      CUBLAS_THREAD_MANAGER<device_t>::synchronize_streams(thread_id_, stream_id_);
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType rhs_device_name>
void Vector<ScalarType, device_name>::set(Vector<ScalarType, rhs_device_name>& rhs) {
  resizeNoCopy(rhs.size_);
  COPY_FROM<rhs_device_name, device_name>::execute(rhs.data_, data_, size_);
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType rhs_device_name>
void Vector<ScalarType, device_name>::set(Vector<ScalarType, rhs_device_name>& rhs,
                                          copy_concurrency_type copy_t) {
  const static DeviceType device_t =
      LIN_ALG::CUBLAS_DEVICE_NAME<rhs_device_name, device_name>::device_t;

  resizeNoCopy(rhs.size());
  assert(thread_id_ > -1 and stream_id_ > -1);

  switch (copy_t) {
    case SYNCHRONOUS:
      COPY_FROM<rhs_device_name, device_name>::execute(rhs.data_, data_, size_);
      break;

    case ASYNCHRONOUS:
      CUBLAS_THREAD_MANAGER<device_t>::synchronize_streams(thread_id_, stream_id_);

      COPY_FROM<rhs_device_name, device_name>::execute(rhs.data_, data_, size_, thread_id_,
                                                       stream_id_);

      CUBLAS_THREAD_MANAGER<device_t>::synchronize_streams(thread_id_, stream_id_);
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

template <typename ScalarType, DeviceType device_name>
void Vector<ScalarType, device_name>::resize(size_t new_size) {
  if (new_size > capacity_) {
    // CUBLAS_THREAD_MANAGER<device_name>::synchronize_streams(thread_id, stream_id);

    int new_capacity = (new_size / 64 + 1) * 64;

    ValueType* new_data = nullptr;

    MEMORY_MANAGEMENT<device_name>::allocate(new_data, new_capacity);
    COPY_FROM<device_name, device_name>::execute(&data_[0], &new_data[0], size_);
    MEMORY_MANAGEMENT<device_name>::deallocate(data_);

    data_ = new_data;
    capacity_ = new_capacity;
    size_ = new_size;

    // CUBLAS_THREAD_MANAGER<device_name>::synchronize_streams(thread_id, stream_id);
  }
  else
    size_ = new_size;
}

template <typename ScalarType, DeviceType device_name>
void Vector<ScalarType, device_name>::resizeNoCopy(size_t new_size) {
  if (new_size > capacity_) {
    // CUBLAS_THREAD_MANAGER<device_name>::synchronize_streams(thread_id, stream_id);

    int new_capacity = (new_size / 64 + 1) * 64;

    MEMORY_MANAGEMENT<device_name>::deallocate(data_);
    MEMORY_MANAGEMENT<device_name>::allocate(data_, new_capacity);

    capacity_ = new_capacity;
    size_ = new_size;

    // CUBLAS_THREAD_MANAGER<device_name>::synchronize_streams(thread_id, stream_id);
  }
  else
    size_ = new_size;
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType rhs_device_name>
ScalarType Vector<ScalarType, device_name>::difference(Vector<ScalarType, rhs_device_name>& rhs) const {
  if (size_ != rhs.size())
    throw std::logic_error(__FUNCTION__);

  Vector<ScalarType, CPU> cp_this(*this);
  Vector<ScalarType, CPU> cp_rhs(rhs);

  ScalarType max_dif = 0;

  for (int i = 0; i < size_; ++i)
    max_dif = std::max(max_dif, std::fabs(cp_this[i] - cp_rhs[i]));

  if (max_dif > 1.e-6)
    throw std::logic_error(__FUNCTION__);

  return max_dif;
}

template <typename ScalarType, DeviceType device_name>
ScalarType Vector<ScalarType, device_name>::difference(std::vector<ScalarType>& rhs) const {
  if (size_ != rhs.size())
    throw std::logic_error(__FUNCTION__);

  Vector<ScalarType, CPU> cp_this(*this);

  ScalarType max_dif = 0;

  for (int i = 0; i < size_; ++i)
    max_dif = std::max(max_dif, std::fabs(cp_this[i] - rhs[i]));

  if (std::fabs(max_dif) > 1.e-6)
    throw std::logic_error(__FUNCTION__);

  return max_dif;
}

template <typename ScalarType, DeviceType device_name>
void Vector<ScalarType, device_name>::print() {
  MEMORY_MANAGEMENT<device_name>::print(data_, size_, capacity_);
}

}  // linalg
}  // dca

#endif  // DCA_LINALG_VECTOR_HPP
