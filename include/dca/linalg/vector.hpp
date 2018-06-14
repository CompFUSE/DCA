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
// This file provides the Vector object for different device types.

#ifndef DCA_LINALG_VECTOR_HPP
#define DCA_LINALG_VECTOR_HPP

#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "dca/linalg/device_type.hpp"
#include "dca/linalg/util/copy.hpp"
#include "dca/linalg/util/memory.hpp"

namespace dca {
namespace linalg {
// dca::linalg::

template <typename ScalarType, DeviceType device_name>
class Vector {
public:
  using ThisType = Vector<ScalarType, device_name>;
  using ValueType = ScalarType;

  Vector();
  Vector(std::string name);

  Vector(size_t size);
  Vector(std::string name, size_t size);

  // Preconditions: capacity >= size.
  Vector(size_t size, size_t capacity);
  Vector(std::string name, size_t size, size_t capacity);

  Vector(const Vector<ScalarType, device_name>& rhs);

  template <DeviceType rhs_device_name>
  Vector(const Vector<ScalarType, rhs_device_name>& rhs);

  ~Vector();

  Vector<ScalarType, device_name>& operator=(const Vector<ScalarType, device_name>& rhs);

  template <DeviceType rhs_device_name>
  Vector<ScalarType, device_name>& operator=(const Vector<ScalarType, rhs_device_name>& rhs);

  Vector<ScalarType, device_name>& operator=(const std::vector<ScalarType>& rhs);

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

  // Asynchronous assignment (copy with stream = getStream(thread_id, stream_id))
  // + synchronization of stream
  // Preconditions: 0 <= thread_id < DCA_MAX_THREADS,
  //                0 <= stream_id < DCA_STREAMS_PER_THREADS.
  void set(const std::vector<ScalarType>& rhs, int thread_id, int stream_id);
  template <DeviceType rhs_device_name>
  void set(const Vector<ScalarType, rhs_device_name>& rhs, int thread_id, int stream_id);

#ifdef DCA_HAVE_CUDA
  template <class Allocator>
  // Asynchronous assignment.
  void setAsync(const std::vector<ScalarType, Allocator>& rhs, cudaStream_t stream);
  template <DeviceType rhs_device>
  void setAsync(const Vector<ScalarType, rhs_device>& rhs, cudaStream_t stream);
#endif  // DCA_HAVE_CUDA

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

  // Releases the memory allocated by *this and sets size and capacity to zero
  void clear();

  // Prints the values of the vector elements.
  template <DeviceType dn = device_name>
  std::enable_if_t<device_name == CPU && dn == CPU, void> print() const;
  template <DeviceType dn = device_name>
  std::enable_if_t<device_name != CPU && dn == device_name, void> print() const;
  // Prints the properties of *this.
  void printFingerprint() const;

private:
  std::string name_;

  size_t size_;
  size_t capacity_;

  ValueType* data_;

  template <typename ScalarType2, DeviceType device_name2>
  friend class dca::linalg::Vector;
};

template <typename ScalarType, DeviceType device_name>
Vector<ScalarType, device_name>::Vector() : Vector("unnamed Vector", 0, 0) {}

template <typename ScalarType, DeviceType device_name>
Vector<ScalarType, device_name>::Vector(std::string name) : Vector(name, 0, 0) {}

template <typename ScalarType, DeviceType device_name>
Vector<ScalarType, device_name>::Vector(size_t size) : Vector("unnamed Vector", size, size) {}

template <typename ScalarType, DeviceType device_name>
Vector<ScalarType, device_name>::Vector(std::string name, size_t size) : Vector(name, size, size) {}

template <typename ScalarType, DeviceType device_name>
Vector<ScalarType, device_name>::Vector(size_t size, size_t capacity)
    : Vector("unnamed Vector", size, capacity) {}

template <typename ScalarType, DeviceType device_name>
Vector<ScalarType, device_name>::Vector(std::string name, size_t size, size_t capacity)
    : name_(name), size_(size), capacity_(capacity), data_(nullptr) {
  assert(capacity_ >= size_);
  util::Memory<device_name>::allocate(data_, capacity_);
  util::Memory<device_name>::setToZero(data_, capacity_);
}

template <typename ScalarType, DeviceType device_name>
Vector<ScalarType, device_name>::Vector(const Vector<ScalarType, device_name>& rhs)
    : name_(rhs.name_), size_(rhs.size_), capacity_(rhs.capacity_), data_(nullptr) {
  assert(capacity_ >= size_);
  util::Memory<device_name>::allocate(data_, capacity_);
  util::memoryCopy(data_, rhs.data_, size_);
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType rhs_device_name>
Vector<ScalarType, device_name>::Vector(const Vector<ScalarType, rhs_device_name>& rhs)
    : name_(rhs.name_), size_(rhs.size_), capacity_(rhs.capacity_), data_(nullptr) {
  assert(capacity_ >= size_);
  util::Memory<device_name>::allocate(data_, capacity_);
  util::memoryCopy(data_, rhs.data_, size_);
}

template <typename ScalarType, DeviceType device_name>
Vector<ScalarType, device_name>::~Vector() {
  util::Memory<device_name>::deallocate(data_);
}

template <typename ScalarType, DeviceType device_name>
Vector<ScalarType, device_name>& Vector<ScalarType, device_name>::operator=(
    const Vector<ScalarType, device_name>& rhs) {
  resizeNoCopy(rhs.size_);

  util::memoryCopy(data_, rhs.data_, size_);
  return *this;
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType rhs_device_name>
Vector<ScalarType, device_name>& Vector<ScalarType, device_name>::operator=(
    const Vector<ScalarType, rhs_device_name>& rhs) {
  resizeNoCopy(rhs.size_);

  util::memoryCopy(data_, rhs.data_, size_);
  return *this;
}

template <typename ScalarType, DeviceType device_name>
Vector<ScalarType, device_name>& Vector<ScalarType, device_name>::operator=(
    const std::vector<ScalarType>& rhs) {
  resizeNoCopy(rhs.size());

  util::memoryCopy(data_, &rhs[0], size_);
  return *this;
}

template <typename ScalarType, DeviceType device_name>
void Vector<ScalarType, device_name>::set(const std::vector<ScalarType>& rhs, int thread_id,
                                          int stream_id) {
  resizeNoCopy(rhs.size());
  util::memoryCopy(data_, &rhs[0], size_, thread_id, stream_id);
}

#ifdef DCA_HAVE_CUDA
template <typename ScalarType, DeviceType device_name>
template <class Allocator>
void Vector<ScalarType, device_name>::setAsync(const std::vector<ScalarType, Allocator>& rhs,
                                               const cudaStream_t stream) {
  resizeNoCopy(rhs.size());
  util::memoryCopyAsync(data_, rhs.data(), size_, stream);
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType rhs_device>
void Vector<ScalarType, device_name>::setAsync(const Vector<ScalarType, rhs_device>& rhs,
                                               const cudaStream_t stream) {
  resizeNoCopy(rhs.size());
  util::memoryCopyAsync(data_, rhs.ptr(), size_, stream);
}
#endif  // DCA_HAVE_CUDA

template <typename ScalarType, DeviceType device_name>
template <DeviceType rhs_device_name>
void Vector<ScalarType, device_name>::set(const Vector<ScalarType, rhs_device_name>& rhs,
                                          int thread_id, int stream_id) {
  resizeNoCopy(rhs.size());
  util::memoryCopy(data_, rhs.data_, size_, thread_id, stream_id);
}

template <typename ScalarType, DeviceType device_name>
void Vector<ScalarType, device_name>::resize(size_t new_size) {
  if (new_size > capacity_) {
    int new_capacity = (new_size / 64 + 1) * 64;

    ValueType* new_data = nullptr;

    util::Memory<device_name>::allocate(new_data, new_capacity);
    util::memoryCopy(new_data, data_, size_);
    util::Memory<device_name>::deallocate(data_);

    data_ = new_data;
    capacity_ = new_capacity;
    size_ = new_size;
  }
  else
    size_ = new_size;
}

template <typename ScalarType, DeviceType device_name>
void Vector<ScalarType, device_name>::resizeNoCopy(size_t new_size) {
  if (new_size > capacity_) {
    int new_capacity = (new_size / 64 + 1) * 64;

    util::Memory<device_name>::deallocate(data_);
    util::Memory<device_name>::allocate(data_, new_capacity);

    capacity_ = new_capacity;
    size_ = new_size;
  }
  else
    size_ = new_size;
}

template <typename ScalarType, DeviceType device_name>
void Vector<ScalarType, device_name>::clear() {
  util::Memory<device_name>::deallocate(data_);
  size_ = capacity_ = 0;
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType dn>
std::enable_if_t<device_name == CPU && dn == CPU, void> Vector<ScalarType, device_name>::print() const {
  printFingerprint();

  std::stringstream ss;
  ss.precision(6);
  ss << std::scientific;

  ss << "\n";
  for (int i = 0; i < size_; i++)
    ss << "\t" << operator[](i);
  ss << "\n";

  std::cout << ss.str() << std::endl;
}

template <typename ScalarType, DeviceType device_name>
template <DeviceType dn>
std::enable_if_t<device_name != CPU && dn == device_name, void> Vector<ScalarType, device_name>::print()
    const {
  Vector<ScalarType, CPU> copy(*this);
  copy.print();
}

template <typename ScalarType, DeviceType device_name>
void Vector<ScalarType, device_name>::printFingerprint() const {
  std::stringstream ss;

  ss << "\n";
  ss << "    name: " << name_ << "\n";
  ss << "    size: " << size_ << "\n";
  ss << "    capacity: " << capacity_ << "\n";
  ss << "    memory-size: " << capacity_ * sizeof(ScalarType) * 1.e-6 << "(Mbytes)\n";

  std::cout << ss.str() << std::endl;
}

}  // linalg
}  // dca

#endif  // DCA_LINALG_VECTOR_HPP
