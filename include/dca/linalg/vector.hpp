// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
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
#include "dca/linalg/util/allocators/allocators.hpp"
#include "dca/linalg/util/copy.hpp"

namespace dca {
namespace linalg {
// dca::linalg::

template <typename ScalarType, DeviceType device_name,
          class Allocator = util::DefaultAllocator<ScalarType, device_name>>
class Vector {
public:
  using ThisType = Vector<ScalarType, device_name, Allocator>;
  using ValueType = ScalarType;

  Vector();
  Vector(const std::string& name);
  Vector(size_t size);
  Vector(const std::string& name, size_t size);
  // Preconditions: capacity >= size.
  Vector(size_t size, size_t capacity);
  Vector(const std::string& name, size_t size, size_t capacity);

  Vector(const ThisType& rhs, const std::string& name = default_name_);

  template <DeviceType device_name2, class Allocator2>
  Vector(const Vector<ScalarType, device_name2, Allocator2>& rhs,
         const std::string& name = default_name_);

  template <DeviceType device_name2, class Allocator2>
  Vector(Vector<ScalarType, device_name2, Allocator2>&& rhs, const std::string& name = default_name_);

  ~Vector();

  ThisType& operator=(const ThisType& rhs);

  template <DeviceType device_name2, class Allocator2>
  ThisType& operator=(const Vector<ScalarType, device_name2, Allocator2>& rhs);

  template <DeviceType device_name2, class Allocator2>
  ThisType& operator=(Vector<ScalarType, device_name2, Allocator2>&& rhs);

  // Returns the i-th element of the vector.
  // Preconditions: 0 <= i < size().first.
  // This method is available only if device_name == CPU.
  template <DeviceType dn = device_name, typename = std::enable_if_t<dn == CPU>>
  ScalarType& operator[](size_t i) {
    assert(i < size_);
    return data_[i];
  }
  template <DeviceType dn = device_name, typename = std::enable_if_t<dn == CPU>>
  const ScalarType& operator[](size_t i) const {
    assert(i < size_);
    return data_[i];
  }

  const std::string& get_name() const {
    return name_;
  }

  void set_name(const std::string& name) {
    name_ = name;
  }

  // Asynchronous assignment (copy with stream = getStream(thread_id, stream_id))
  // + synchronization of stream.
  template <class Container>
  void set(const Container& rhs, int thread_id, int stream_id);

  // Synchronous copy. Container must define a data() and size() methods.
  // Precondition: rhs has the same size as this Vector.
  template <class Container>
  void copyTo(Container& rhs) const;

#ifdef DCA_HAVE_CUDA
  // Asynchronous assignment.
  template <class Container>
  void setAsync(const Container& rhs, cudaStream_t stream);

  void setToZeroAsync(cudaStream_t stream);
#endif  // DCA_HAVE_CUDA
  template <class Container>
  void setAsync(const Container& rhs, int thred_id, int stream_id = 0);

  // Returns the pointer to the 0-th element of the vector.
  ValueType* ptr() {
    return data_;
  }
  ValueType* data() {
    return data_;
  }
  const ValueType* ptr() const {
    return data_;
  }
  const ValueType* data() const {
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
  void print() const;
  // Prints the properties of *this.
  void printFingerprint() const;

protected:
  Allocator allocator_;

private:
  std::string name_;

  size_t size_;
  size_t capacity_;

  ValueType* data_;

  template <typename ScalarType2, DeviceType device_name2, class Allocator2>
  friend class dca::linalg::Vector;

  static const std::string default_name_;
};
template <typename ScalarType, DeviceType device_name, class Allocator>
const std::string Vector<ScalarType, device_name, Allocator>::default_name_("no-name");

template <typename ScalarType, DeviceType device_name, class Allocator>
Vector<ScalarType, device_name, Allocator>::Vector() : Vector(default_name_, 0, 0) {}

template <typename ScalarType, DeviceType device_name, class Allocator>
Vector<ScalarType, device_name, Allocator>::Vector(const std::string& name) : Vector(name, 0, 0) {}

template <typename ScalarType, DeviceType device_name, class Allocator>
Vector<ScalarType, device_name, Allocator>::Vector(size_t size)
    : Vector(default_name_, size, size) {}

template <typename ScalarType, DeviceType device_name, class Allocator>
Vector<ScalarType, device_name, Allocator>::Vector(const std::string& name, size_t size)
    : Vector(name, size, size) {}

template <typename ScalarType, DeviceType device_name, class Allocator>
Vector<ScalarType, device_name, Allocator>::Vector(size_t size, size_t capacity)
    : Vector(default_name_, size, capacity) {}

template <typename ScalarType, DeviceType device_name, class Allocator>
Vector<ScalarType, device_name, Allocator>::Vector(const std::string& name, size_t size,
                                                   size_t capacity)
    : name_(name), size_(size), capacity_(capacity), data_(nullptr) {
  assert(capacity_ >= size_);
  data_ = allocator_.allocate(capacity_);
}

template <typename ScalarType, DeviceType device_name, class Allocator>
Vector<ScalarType, device_name, Allocator>::Vector(const ThisType& rhs, const std::string& name)
    : Vector(name) {
  *this = rhs;
}

template <typename ScalarType, DeviceType device_name, class Allocator>
template <DeviceType device_name2, class Allocator2>
Vector<ScalarType, device_name, Allocator>::Vector(
    const Vector<ScalarType, device_name2, Allocator2>& rhs, const std::string& name)
    : Vector(name) {
  *this = rhs;
}

template <typename ScalarType, DeviceType device_name, class Allocator>
template <DeviceType device_name2, class Allocator2>
Vector<ScalarType, device_name, Allocator>::Vector(Vector<ScalarType, device_name2, Allocator2>&& rhs,
                                                   const std::string& name)
    : Vector(name) {
  *this = std::move(rhs);
}

template <typename ScalarType, DeviceType device_name, class Allocator>
Vector<ScalarType, device_name, Allocator>::~Vector() {
  allocator_.deallocate(data_);
}

template <typename ScalarType, DeviceType device_name, class Allocator>
Vector<ScalarType, device_name, Allocator>& Vector<ScalarType, device_name, Allocator>::operator=(
    const ThisType& rhs) {
  resizeNoCopy(rhs.size());
  if (device_name == CPU)
    util::memoryCopyCpu(data_, rhs.data_, size_);
  else
    util::memoryCopy(data_, rhs.data_, size_);

  return *this;
}

template <typename ScalarType, DeviceType device_name, class Allocator>
template <DeviceType device_name2, class Allocator2>
Vector<ScalarType, device_name, Allocator>& Vector<ScalarType, device_name, Allocator>::operator=(
    const Vector<ScalarType, device_name2, Allocator2>& rhs) {
  resizeNoCopy(rhs.size());
  if (device_name == CPU && device_name2 == CPU)
    util::memoryCopyCpu(data_, rhs.data_, size_);
  else
    util::memoryCopy(data_, rhs.data_, size_);

  return *this;
}

template <typename ScalarType, DeviceType device_name, class Allocator>
template <DeviceType device_name2, class Allocator2>
Vector<ScalarType, device_name, Allocator>& Vector<ScalarType, device_name, Allocator>::operator=(
    Vector<ScalarType, device_name2, Allocator2>&& rhs) {
  std::swap(data_, rhs.data_);
  std::swap(size_, rhs.size_);
  std::swap(capacity_, rhs.capacity_);

  return *this;
}

template <typename ScalarType, DeviceType device_name, class Allocator>
template <class Container>
void Vector<ScalarType, device_name, Allocator>::set(const Container& rhs, int thread_id,
                                                     int stream_id) {
  resizeNoCopy(rhs.size());
  util::memoryCopy(data_, rhs.data(), size_, thread_id, stream_id);
}

template <typename ScalarType, DeviceType device_name, class Allocator>
template <class Container>
void Vector<ScalarType, device_name, Allocator>::copyTo(Container& rhs) const {
  if (rhs.size() != size())
    throw(std::logic_error("The size of the destination container is different."));
  util::memoryCopy(rhs.data(), data_, size_);
}

#ifdef DCA_HAVE_CUDA
template <typename ScalarType, DeviceType device_name, class Allocator>
template <class Container>
void Vector<ScalarType, device_name, Allocator>::setAsync(const Container& rhs,
                                                          const cudaStream_t stream) {
  resizeNoCopy(rhs.size());
  util::memoryCopyAsync(data_, rhs.data(), size_, stream);
}

template <typename ScalarType, DeviceType device_name, class Allocator>
void Vector<ScalarType, device_name, Allocator>::setToZeroAsync(cudaStream_t stream) {
  cudaMemsetAsync(data_, 0, size_ * sizeof(ScalarType), stream);
}
#endif  // DCA_HAVE_CUDA

template <typename ScalarType, DeviceType device_name, class Allocator>
template <class Container>
void Vector<ScalarType, device_name, Allocator>::setAsync(const Container& rhs, const int thread_id,
                                                          const int stream_id) {
#ifdef DCA_HAVE_CUDA
  setAsync(rhs, util::getStream(thread_id, stream_id));
#else
  set(rhs, thread_id, stream_id);
#endif  // DCA_HAVE_CUDA
}

template <typename ScalarType, DeviceType device_name, class Allocator>
void Vector<ScalarType, device_name, Allocator>::resize(size_t new_size) {
  if (new_size > capacity_) {
    int new_capacity = (new_size / 64 + 1) * 64;

    ValueType* new_data = allocator_.allocate(new_capacity);
    util::memoryCopy(new_data, data_, size_);
    allocator_.deallocate(data_);

    data_ = new_data;
    capacity_ = new_capacity;
    size_ = new_size;
  }
  else
    size_ = new_size;
}

template <typename ScalarType, DeviceType device_name, class Allocator>
void Vector<ScalarType, device_name, Allocator>::resizeNoCopy(size_t new_size) {
  if (new_size > capacity_) {
    int new_capacity = (new_size / 64 + 1) * 64;

    allocator_.deallocate(data_);
    data_ = allocator_.allocate(new_capacity);

    capacity_ = new_capacity;
    size_ = new_size;
  }
  else
    size_ = new_size;
}

template <typename ScalarType, DeviceType device_name, class Allocator>
void Vector<ScalarType, device_name, Allocator>::clear() {
  allocator_.deallocate(data_);
  size_ = capacity_ = 0;
}

template <typename ScalarType, DeviceType device_name, class Allocator>
void Vector<ScalarType, device_name, Allocator>::print() const {
  if (device_name == GPU) {
    Vector<ScalarType, CPU> copy(*this);
    return copy.print();
  }

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

template <typename ScalarType, DeviceType device_name, class Allocator>
void Vector<ScalarType, device_name, Allocator>::printFingerprint() const {
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
