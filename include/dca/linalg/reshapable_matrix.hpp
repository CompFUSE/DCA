// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides a matrix with a more efficient reshaping between different rectangular shapes
//  of similar total size.

#pragma once

#include <cassert>
#include <cmath>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include "dca/linalg/util/allocators/allocators.hpp"
#include "dca/linalg/device_type.hpp"
#include "dca/linalg/util/copy.hpp"
#include "dca/linalg/util/stream_functions.hpp"

namespace dca {
namespace linalg {
// dca::linalg::

template <typename ScalarType, DeviceType device_name,
          class Allocator = util::DefaultAllocator<ScalarType, device_name>>
class ReshapableMatrix : public Allocator {
public:
  using ThisType = ReshapableMatrix<ScalarType, device_name, Allocator>;
  using ValueType = ScalarType;

  // Default contructor creates a matrix of zero size and capacity.
  ReshapableMatrix() = default;
  // Initializes a square size x size matrix.
  ReshapableMatrix(int size);
  // Initializes a square size.first x size.second matrix.
  ReshapableMatrix(std::pair<int, int> size);

  // Contructs a matrix with size rhs.size() and a copy of the elements of rhs.
  ReshapableMatrix(const ThisType& rhs);
  template <DeviceType rhs_device_name, class AllocatorRhs>
  ReshapableMatrix(const ReshapableMatrix<ScalarType, rhs_device_name, AllocatorRhs>& rhs);

  // Constructs a matrix with size rhs.size(). The elements of rhs are moved.
  ReshapableMatrix(ThisType&& rhs);

  // Resize the matrix to rhs.size() and copies the elements.
  ReshapableMatrix& operator=(const ThisType& rhs);
  template <DeviceType rhs_device_name, class AllocatorRhs>
  ReshapableMatrix& operator=(const ReshapableMatrix<ScalarType, rhs_device_name, AllocatorRhs>& rhs);

  // Moves the elements of rhs into this matrix.
  ReshapableMatrix& operator=(ThisType&& rhs);

  ~ReshapableMatrix();

  // Returns true if this is equal to other, false otherwise.
  // Two matrices are equal, if they have the same size and contain the same elements. Name and
  // capacity are ignored.
  // Special case: two matrices without elements are equal.
  bool operator==(const ReshapableMatrix<ScalarType, device_name, Allocator>& other) const;

  // Returns true if this is not equal to other, false otherwise.
  // See description of operator== for the definition of equality.
  bool operator!=(const ReshapableMatrix<ScalarType, device_name, Allocator>& other) const;

  // Returns the (i,j)-th element of the matrix.
  // Preconditions: 0 <= i < size().first, 0 <= j < size().second.
  // This method is available only if device_name == CPU.
  template <DeviceType dn = device_name, typename = std::enable_if_t<dn == CPU>>
  ScalarType& operator()(int i, int j) {
    assert(i >= 0 && i < size_.first);
    assert(j >= 0 && j < size_.second);
    return data_[i + j * leadingDimension()];
  }
  template <DeviceType dn = device_name, typename = std::enable_if_t<dn == CPU>>
  const ScalarType& operator()(int i, int j) const {
    assert(i >= 0 && i < size_.first);
    assert(j >= 0 && j < size_.second);
    return data_[i + j * leadingDimension()];
  }

  // Returns the pointer to the (0,0)-th element.
  ValueType* ptr() {
    return data_;
  }
  const ValueType* ptr() const {
    return data_;
  }

  // Returns the pointer to the (i,j)-th element i < size().first and 0 < j < size().second, or
  // a pointer past the end of the range if i == size().first or j == size().second.
  // Preconditions: 0 <= i <= size().first, 0 <= j <= size().second.
  ValueType* ptr(int i, int j) {
    assert(i >= 0 && i <= size_.first);
    assert(j >= 0 && j <= size_.second);
    return data_ + i + j * leadingDimension();
  }
  const ValueType* ptr(int i, int j) const {
    assert(i >= 0 && i <= size_.first);
    assert(j >= 0 && j <= size_.second);
    return data_ + i + j * leadingDimension();
  }

  const std::pair<int, int> size() const {
    return size_;
  }
  std::size_t capacity() const {
    return capacity_;
  }
  int nrRows() const {
    return size_.first;
  }
  int nrCols() const {
    return size_.second;
  }
  int leadingDimension() const {
    return size_.first;
  }

  // Resizes *this to a (new_size.first * new_size.second) matrix.
  // The previous elements are not copied, therefore all the elements
  // may have any value after the call to this method.
  // Returns: true if reallocation took place.
  bool resizeNoCopy(std::pair<int, int> new_size);
  // Resizes *this to a (new_size * new_size) matrix. See previous method for details.
  bool resizeNoCopy(int new_size) {
    return resizeNoCopy(std::make_pair(new_size, new_size));
  }
  // Resizes *this to a new size of rhs matrix and perform element-wise copy from rhs
  void copyFrom(const ThisType& rhs);

  // Reserves the space for at least (new_size.first * new_size.second) elements without changing
  // the matrix size. The value of the matrix elements is undefined after calling this method.
  // Returns: true if reallocation took place.
  bool reserveNoCopy(std::size_t new_size);

  void swap(ReshapableMatrix<ScalarType, device_name, Allocator>& other);

  // Releases the memory allocated by *this and sets size and capacity to zero.
  void clear();

#ifdef DCA_HAVE_CUDA
  // Asynchronous assignment.
  template <DeviceType rhs_device_name>
  void setAsync(const ReshapableMatrix<ScalarType, rhs_device_name>& rhs, cudaStream_t stream);

  // Asynchronous assignment (copy with stream = getStream(thread_id, stream_id))
  template <DeviceType rhs_device_name>
  void setAsync(const ReshapableMatrix<ScalarType, rhs_device_name>& rhs, int thread_id,
                int stream_id);

  void setToZero(cudaStream_t stream);
#else
  // Synchronous assignment fallback for SetAsync.
  template <DeviceType rhs_device_name>
  void setAsync(const ReshapableMatrix<ScalarType, rhs_device_name>& rhs, int /*thread_id*/,
                int /*stream_id*/);

#endif  // DCA_HAVE_CUDA

  // Returns the allocated device memory in bytes.
  std::size_t deviceFingerprint() const;

private:
  static std::size_t nextCapacity(std::size_t size);
  inline static size_t nrElements(std::pair<int, int> size) {
    return static_cast<size_t>(size.first) * static_cast<size_t>(size.second);
  }

  std::pair<int, int> size_ = std::make_pair(0, 0);
  std::size_t capacity_ = 0;

  ValueType* data_ = nullptr;

  template <class ScalarType2, DeviceType device_name2, class Allocator2>
  friend class dca::linalg::ReshapableMatrix;
};

template <typename ScalarType, DeviceType device_name, class Allocator>
ReshapableMatrix<ScalarType, device_name, Allocator>::ReshapableMatrix(int size)
    : ReshapableMatrix(std::make_pair(size, size)) {}

template <typename ScalarType, DeviceType device_name, class Allocator>
ReshapableMatrix<ScalarType, device_name, Allocator>::ReshapableMatrix(std::pair<int, int> size)
    : size_(size), capacity_(nextCapacity(nrElements(size))) {
  assert(size_.first >= 0 && size_.second >= 0);
  assert(capacity_ >= nrElements(size_));

  data_ = Allocator::allocate(capacity_);
}

template <typename ScalarType, DeviceType device_name, class Allocator>
ReshapableMatrix<ScalarType, device_name, Allocator>::ReshapableMatrix(const ThisType& rhs) {
  *this = rhs;
}

template <typename ScalarType, DeviceType device_name, class Allocator>
template <DeviceType rhs_device_name, class AllocatorRhs>
ReshapableMatrix<ScalarType, device_name, Allocator>::ReshapableMatrix(
    const ReshapableMatrix<ScalarType, rhs_device_name, AllocatorRhs>& rhs) {
  *this = rhs;
}

template <typename ScalarType, DeviceType device_name, class Allocator>
ReshapableMatrix<ScalarType, device_name, Allocator>::ReshapableMatrix(
    ReshapableMatrix<ScalarType, device_name, Allocator>&& rhs)
    : ReshapableMatrix<ScalarType, device_name, Allocator>() {
  swap(rhs);
}

template <typename ScalarType, DeviceType device_name, class Allocator>
ReshapableMatrix<ScalarType, device_name, Allocator>& ReshapableMatrix<
    ScalarType, device_name, Allocator>::operator=(const ThisType& rhs) {
  size_ = rhs.size_;
  capacity_ = rhs.capacity_;

  Allocator::deallocate(data_);
  data_ = Allocator::allocate(capacity_);
  util::memoryCopy(data_, leadingDimension(), rhs.data_, rhs.leadingDimension(), size_);
  return *this;
}

template <typename ScalarType, DeviceType device_name, class Allocator>
template <DeviceType rhs_device_name, class AllocatorRhs>
ReshapableMatrix<ScalarType, device_name, Allocator>& ReshapableMatrix<
    ScalarType, device_name,
    Allocator>::operator=(const ReshapableMatrix<ScalarType, rhs_device_name, AllocatorRhs>& rhs) {
  size_ = rhs.size_;
  capacity_ = rhs.capacity_;

  Allocator::deallocate(data_);
  data_ = Allocator::allocate(capacity_);
  util::memoryCopy(data_, leadingDimension(), rhs.data_, rhs.leadingDimension(), size_);
  return *this;
}

template <typename ScalarType, DeviceType device_name, class Allocator>
ReshapableMatrix<ScalarType, device_name, Allocator>& ReshapableMatrix<
    ScalarType, device_name, Allocator>::operator=(ThisType&& rhs) {
  swap(rhs);
  return *this;
}

template <typename ScalarType, DeviceType device_name, class Allocator>
ReshapableMatrix<ScalarType, device_name, Allocator>::~ReshapableMatrix() {
  Allocator::deallocate(data_);
}

template <typename ScalarType, DeviceType device_name, class Allocator>
bool ReshapableMatrix<ScalarType, device_name, Allocator>::operator==(
    const ReshapableMatrix<ScalarType, device_name, Allocator>& other) const {
  if (device_name == GPU)
    return ReshapableMatrix<ScalarType, CPU>(*this) == ReshapableMatrix<ScalarType, CPU>(other);

  if (size() != other.size())
    return nrRows() * nrCols() == 0 and other.nrRows() * other.nrCols() == 0;

  for (int j = 0; j < nrCols(); ++j)
    for (int i = 0; i < nrRows(); ++i)
      if ((*this)(i, j) != other(i, j))
        return false;

  return true;
}

template <typename ScalarType, DeviceType device_name, class Allocator>
bool ReshapableMatrix<ScalarType, device_name, Allocator>::operator!=(
    const ReshapableMatrix<ScalarType, device_name, Allocator>& other) const {
  return !(*this == other);
}

template <typename ScalarType, DeviceType device_name, class Allocator>
bool ReshapableMatrix<ScalarType, device_name, Allocator>::resizeNoCopy(
    const std::pair<int, int> new_size) {
  const bool realloc = reserveNoCopy(nrElements(new_size));
  size_ = new_size;
  assert(capacity_ >= nrElements(size_));
  return realloc;
}

template <typename ScalarType, DeviceType device_name, class Allocator>
bool ReshapableMatrix<ScalarType, device_name, Allocator>::reserveNoCopy(std::size_t new_size) {
  if (new_size > capacity_) {
    Allocator::deallocate(data_);
    capacity_ = nextCapacity(new_size);
    data_ = Allocator::allocate(capacity_);
    return true;
  }
  return false;
}

template <typename ScalarType, DeviceType device_name, class Allocator>
void ReshapableMatrix<ScalarType, device_name, Allocator>::copyFrom(const ThisType& rhs) {
  auto new_size = rhs.size();
  auto new_nrElements = nrElements(new_size);

  if (new_nrElements > capacity_) {
     Allocator::deallocate(data_);
     capacity_ = nextCapacity(new_nrElements);
     data_ = Allocator::allocate(capacity_);
  }

  util::memoryCopy(data_, leadingDimension(), rhs.data_, rhs.leadingDimension(), size_);
}

template <typename ScalarType, DeviceType device_name, class Allocator>
void ReshapableMatrix<ScalarType, device_name, Allocator>::swap(
    ReshapableMatrix<ScalarType, device_name, Allocator>& other) {
  std::swap(size_, other.size_);
  std::swap(capacity_, other.capacity_);
  std::swap(data_, other.data_);
}

template <typename ScalarType, DeviceType device_name, class Allocator>
void ReshapableMatrix<ScalarType, device_name, Allocator>::clear() {
  Allocator::deallocate(data_);
  size_ = std::make_pair(0, 0);
  capacity_ = 0;
}

#ifdef DCA_HAVE_CUDA

template <typename ScalarType, DeviceType device_name, class Allocator>
template <DeviceType rhs_device_name>
void ReshapableMatrix<ScalarType, device_name, Allocator>::setAsync(
    const ReshapableMatrix<ScalarType, rhs_device_name>& rhs, const cudaStream_t stream) {
  resizeNoCopy(rhs.size_);
  util::memoryCopyAsync(data_, leadingDimension(), rhs.data_, rhs.leadingDimension(), size_, stream);
}

template <typename ScalarType, DeviceType device_name, class Allocator>
template <DeviceType rhs_device_name>
void ReshapableMatrix<ScalarType, device_name, Allocator>::setAsync(
    const ReshapableMatrix<ScalarType, rhs_device_name>& rhs, const int thread_id,
    const int stream_id) {
  setAsync(rhs, util::getStream(thread_id, stream_id));
}

template <typename ScalarType, DeviceType device_name, class Allocator>
void ReshapableMatrix<ScalarType, device_name, Allocator>::setToZero(cudaStream_t stream) {
  cudaMemsetAsync(data_, 0, leadingDimension() * nrCols() * sizeof(ScalarType), stream);
}

#else  // DCA_HAVE_CUDA

template <typename ScalarType, DeviceType device_name, class Allocator>
template <DeviceType rhs_device_name>
void ReshapableMatrix<ScalarType, device_name, Allocator>::setAsync(
    const ReshapableMatrix<ScalarType, rhs_device_name>& rhs, int /*thread_id*/, int /*stream_id*/) {
  *this = rhs;
}

#endif  // DCA_HAVE_CUDA

template <typename ScalarType, DeviceType device_name, class Allocator>
std::size_t ReshapableMatrix<ScalarType, device_name, Allocator>::nextCapacity(const std::size_t size) {
  assert(size >= 0);
  constexpr std::size_t block_size = 512;

  auto next_power_of_two = [](std::size_t x) {
    if (!x)
      return std::size_t(0);
    std::size_t result = 1;
    while (result < x) {
      result <<= 1;
    }
    return result;
  };

  return size <= block_size ? next_power_of_two(size)
                            : (size + block_size - 1) / block_size * block_size;
}

template <typename ScalarType, DeviceType device_name, class Allocator>
std::size_t ReshapableMatrix<ScalarType, device_name, Allocator>::deviceFingerprint() const {
  if (device_name == GPU)
    return capacity_ * sizeof(ScalarType);
  else
    return 0;
}

}  // namespace linalg
}  // namespace dca
