// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//         Peter W. Doak (doakpw@ornl.gov
//
/// \file provides the Matrix object for different device types and allocators.

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

#include "dca/linalg/vector.hpp"
#include "dca/linalg/util/allocators/allocators.hpp"
#include "dca/linalg/device_type.hpp"
#include "dca/linalg/util/copy.hpp"
#include "dca/linalg/util/memory.hpp"
#include "dca/linalg/util/stream_functions.hpp"
#include "dca/util/type_help.hpp"

namespace dca {
namespace linalg {
// dca::linalg::

/** Matrix class for interfacing with Blas, Cublas, Rocblas
 *  its row major i.e, row is fast.
 */
template <typename ScalarType, DeviceType device_name, class ALLOC = util::DefaultAllocator<ScalarType, device_name>>
class Matrix : public ALLOC {
public:
  using ThisType = Matrix<ScalarType, device_name>;
  using ValueType = ScalarType;
  using Allocator = ALLOC;
  constexpr static DeviceType device = device_name;

  Matrix(const std::string& name = default_name_);

  Matrix(int size);
  Matrix(const std::string& name, int size);

  // Preconditions: capacity >= size.
  Matrix(int size, int capacity);
  Matrix(const std::string& name, int size, int capacity);

  Matrix(std::pair<int, int> size);
  Matrix(const std::string& name, std::pair<int, int> size);

  // Preconditions: capacity.first >= size.first, capacity.second >= size.second.
  Matrix(std::pair<int, int> size, std::pair<int, int> capacity);
  Matrix(const std::string& name, std::pair<int, int> size, std::pair<int, int> capacity);

  // Copy and move constructor:
  // Constructs a matrix with name name, size rhs.size() and a copy of the elements of rhs.
  Matrix(const Matrix<ScalarType, device_name, ALLOC>& rhs, const std::string& name = default_name_);
  // Constructs a matrix with name name, size rhs.size(). The elements of rhs are moved.
  // Postcondition: rhs is a (0 x 0) matrix.
  Matrix(Matrix<ScalarType, device_name, ALLOC>&& rhs, const std::string& = default_name_);

  // Contructs a matrix with name name, size rhs.size() and a copy of the elements of rhs, where rhs
  // elements are stored on a different device.
  template <DeviceType rhs_device_name, class rhs_ALLOC>
  Matrix(const Matrix<ScalarType, rhs_device_name, rhs_ALLOC>& rhs, const std::string& = default_name_);

  // Contructs a matrix with name name, size rhs.size() and a copy of the elements of rhs, where rhs
  // elements are stored on a different device.
  template <typename Scalar2, DeviceType rhs_device_name, class rhs_ALLOC>
  Matrix(const Matrix<Scalar2, rhs_device_name, rhs_ALLOC>& rhs, const std::string& = default_name_);

  ~Matrix();

  // Assignment operators:
  // Resizes the matrix to rhs.size() and copy the elements of rhs.
  // Postcondition: The name of the matrix is unchanged.
  Matrix<ScalarType, device_name, ALLOC>& operator=(const Matrix<ScalarType, device_name, ALLOC>& rhs);
  // Resizes the matrix to rhs.size() and move the elements of rhs.
  // Postcondition: The name of the matrix is unchanged; rhs is a (0 x 0) matrix.
  Matrix<ScalarType, device_name, ALLOC>& operator=(Matrix<ScalarType, device_name, ALLOC>&& rhs);

  // Resizes the matrix to rhs.size() and copy the elements, stored on a different device, of rhs.
  // Postcondition: The name of the matrix is unchanged.
  template <DeviceType rhs_device_name, class rhs_ALLOC>
  Matrix<ScalarType, device_name, ALLOC>& operator=(const Matrix<ScalarType, rhs_device_name, rhs_ALLOC>& rhs);

  template <typename ScalarRhs, DeviceType rhs_device_name, class rhs_ALLOC>
  Matrix<ScalarType, device_name, ALLOC>& operator=(const Matrix<ScalarRhs, rhs_device_name, rhs_ALLOC>& rhs);

  // Returns true if this is equal to other, false otherwise.
  // Two matrices are equal, if they have the same size and contain the same elements. Name and
  // capacity are ignored.
  // Special case: two matrices without elements are equal.
  bool operator==(const Matrix<ScalarType, device_name, ALLOC>& other) const;

  // Returns true if this is not equal to other, false otherwise.
  // See description of operator== for the definition of equality.
  bool operator!=(const Matrix<ScalarType, device_name, ALLOC>& other) const;

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

  bool is_square() const {
    return (size_.first == size_.second);
  }

  const std::pair<int, int> size() const {
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

  int getActualSize() {
    return nrElements(capacity_);
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

  // Releases the memory allocated by *this and sets size and capacity to zero.
  void clear();

  // Swaps the contents of the matrix except the name with those of rhs.
  void swap(Matrix<ScalarType, device_name, ALLOC>& rhs);
  // Swaps the contents of the matrix, included the name, with those of rhs.
  void swapWithName(Matrix<ScalarType, device_name, ALLOC>& rhs);

  // Asynchronous assignment (copy with stream = getStream(thread_id, stream_id))
  // + synchronization of stream
  template <DeviceType rhs_device_name, class rhs_ALLOC>
  void set(const Matrix<ScalarType, rhs_device_name, rhs_ALLOC>& rhs, int thread_id, int stream_id);

  template <DeviceType rhs_device_name, class rhs_ALLOC>
  void set(const Matrix<ScalarType, rhs_device_name, rhs_ALLOC>& rhs, const util::GpuStream& stream);

  // Asynchronous assignment.
  template <DeviceType rhs_device_name, class rhs_ALLOC>
  void setAsync(const Matrix<ScalarType, rhs_device_name, rhs_ALLOC>& rhs, const util::GpuStream& stream);

  // Asynchronous assignment (copy with stream = getStream(thread_id, stream_id))
  template <DeviceType rhs_device_name, class rhs_ALLOC>
  void setAsync(const Matrix<ScalarType, rhs_device_name, rhs_ALLOC>& rhs, int thread_id, int stream_id);

  void setToZero(const util::GpuStream& stream);

  // Prints the values of the matrix elements.
  void print() const;
  // Prints the properties of *this.
  void printFingerprint() const;
  // Returns the allocated device memory in bytes.
  std::size_t deviceFingerprint() const;

  std::string toStr() const;
private:
  static std::pair<int, int> capacityMultipleOfBlockSize(std::pair<int, int> size);
  inline static size_t nrElements(std::pair<int, int> size) {
    return static_cast<size_t>(size.first) * static_cast<size_t>(size.second);
  }
  static constexpr int block_size_ = 32;
  static const std::string default_name_;

  std::string name_;

  std::pair<int, int> size_;
  std::pair<int, int> capacity_;

  ValueType* data_ = nullptr;

  template <class ScalarType2, DeviceType device_name2, class ALLOC2>
  friend class dca::linalg::Matrix;
};

template <typename ScalarType, DeviceType device_name, class ALLOC>
const std::string Matrix<ScalarType, device_name,  ALLOC>::default_name_ = "no-name";

template <typename ScalarType, DeviceType device_name, class ALLOC>
Matrix<ScalarType, device_name,  ALLOC>::Matrix(const std::string& name) : Matrix(name, 0) {}

template <typename ScalarType, DeviceType device_name, class ALLOC>
Matrix<ScalarType, device_name,  ALLOC>::Matrix(int size) : Matrix(std::make_pair(size, size)) {}

template <typename ScalarType, DeviceType device_name, class ALLOC>
Matrix<ScalarType, device_name,  ALLOC>::Matrix(const std::string& name, int size)
    : Matrix(name, std::make_pair(size, size)) {}

template <typename ScalarType, DeviceType device_name, class ALLOC>
Matrix<ScalarType, device_name,  ALLOC>::Matrix(int size, int capacity)
    : Matrix(std::make_pair(size, size), std::make_pair(capacity, capacity)) {}

template <typename ScalarType, DeviceType device_name, class ALLOC>
Matrix<ScalarType, device_name,  ALLOC>::Matrix(const std::string& name, int size, int capacity)
    : Matrix(name, std::make_pair(size, size), std::make_pair(capacity, capacity)) {}

template <typename ScalarType, DeviceType device_name, class ALLOC>
Matrix<ScalarType, device_name,  ALLOC>::Matrix(std::pair<int, int> size) : Matrix(size, size) {}

template <typename ScalarType, DeviceType device_name, class ALLOC>
Matrix<ScalarType, device_name,  ALLOC>::Matrix(const std::string& name, std::pair<int, int> size)
    : Matrix(name, size, size) {}

template <typename ScalarType, DeviceType device_name, class ALLOC>
Matrix<ScalarType, device_name,  ALLOC>::Matrix(std::pair<int, int> size, std::pair<int, int> capacity)
    : Matrix(default_name_, size, capacity) {}

template <typename ScalarType, DeviceType device_name, class ALLOC>
template <DeviceType rhs_device_name, class rhs_ALLOC>
Matrix<ScalarType, device_name, ALLOC>::Matrix(const Matrix<ScalarType, rhs_device_name, rhs_ALLOC>& rhs,
                                        const std::string& name)
    : name_(name), size_(rhs.size_), capacity_(rhs.capacity_) {
  data_ = Allocator::allocate(nrElements(capacity_));
  util::memoryCopy(data_, leadingDimension(), rhs.data_, rhs.leadingDimension(), size_);
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
template <typename ScalarRhs, DeviceType rhs_device_name, class rhs_ALLOC>
Matrix<ScalarType, device_name,  ALLOC>::Matrix(const Matrix<ScalarRhs, rhs_device_name, rhs_ALLOC>& rhs,
                                        const std::string& name)
    : name_(name), size_(rhs.size_), capacity_(rhs.capacity_) {
  if (sizeof(ScalarType) != sizeof(ScalarRhs))
    throw std::runtime_error("conversion of both type and location of Matrix not currently possible!");
  data_ = ALLOC::allocate(nrElements(capacity_));
  util::memoryCopy(data_, leadingDimension(), rhs.data_, rhs.leadingDimension(), size_);
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
Matrix<ScalarType, device_name, ALLOC>::Matrix(const std::string& name, std::pair<int, int> size,
                                        std::pair<int, int> capacity)
    : name_(name), size_(size), capacity_(capacityMultipleOfBlockSize(capacity)) {
  assert(size_.first >= 0 && size_.second >= 0);
  assert(capacity.first >= 0 && capacity.second >= 0);
  assert(capacity.first >= size_.first && capacity.second >= size_.second);
  assert(capacity_.first >= capacity.first && capacity_.second >= capacity.second);

  data_ = ALLOC::allocate(nrElements(capacity_));
  util::Memory<device_name>::setToZero(data_, nrElements(capacity_));
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
Matrix<ScalarType, device_name,  ALLOC>::Matrix(const Matrix<ScalarType, device_name,  ALLOC>& rhs,
                                        const std::string& name)
    : name_(name) {
  *this = rhs;
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
Matrix<ScalarType, device_name,  ALLOC>::Matrix(Matrix<ScalarType, device_name,  ALLOC>&& rhs, const std::string& name)
    : name_(name), size_(rhs.size_), capacity_(rhs.capacity_), data_(rhs.data_) {
  rhs.capacity_ = std::make_pair(0, 0);
  rhs.size_ = std::make_pair(0, 0);
  rhs.data_ = nullptr;
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
Matrix<ScalarType, device_name,  ALLOC>::~Matrix() {
  Allocator::deallocate(data_);
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
void Matrix<ScalarType, device_name,  ALLOC>::resize(std::pair<int, int> new_size) {
  if (new_size.first == 0 || new_size.second ==0) {
    size_ = new_size;
    return;
  } else if (new_size.first > capacity_.first || new_size.second > capacity_.second) {
    std::pair<int, int> new_capacity = capacityMultipleOfBlockSize(new_size);

    ValueType* new_data = nullptr;
    new_data = Allocator::allocate(nrElements(new_capacity));
    // hip memorycpy2D routines don't tolerate leadingDimension = 0
    const std::pair<int, int> copy_size(std::min(new_size.first, size_.first),
                                        std::min(new_size.second, size_.second));
    util::memoryCopy(new_data, new_capacity.first, data_, leadingDimension(), copy_size);
    Allocator::deallocate(data_);
    data_ = new_data;
    capacity_ = new_capacity;
    size_ = new_size;
  }
  else {
    size_ = new_size;
  }
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
Matrix<ScalarType, device_name,  ALLOC>& Matrix<ScalarType, device_name, ALLOC>::operator=(
    const Matrix<ScalarType, device_name,  ALLOC>& rhs) {
  resizeNoCopy(rhs.size_);
  if (device_name == CPU)
    util::memoryCopyCpu(data_, leadingDimension(), rhs.data_, rhs.leadingDimension(), size_);
  else
    util::memoryCopy(data_, leadingDimension(), rhs.data_, rhs.leadingDimension(), size_);
  return *this;
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
Matrix<ScalarType, device_name,  ALLOC>& Matrix<ScalarType, device_name,  ALLOC>::operator=(
    Matrix<ScalarType, device_name,  ALLOC>&& rhs) {
  swap(rhs);
  return *this;
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
template <DeviceType rhs_device_name, class rhs_ALLOC>
Matrix<ScalarType, device_name,  ALLOC>& Matrix<ScalarType, device_name,  ALLOC>::operator=(
    const Matrix<ScalarType, rhs_device_name, rhs_ALLOC>& rhs) {
  resizeNoCopy(rhs.size_);
  util::memoryCopy(data_, leadingDimension(), rhs.data_, rhs.leadingDimension(), size_);
  return *this;
}

#ifdef DCA_HAVE_GPU
template <typename ScalarType, DeviceType device_name, class ALLOC>
template <typename ScalarRhs, DeviceType rhs_device_name, class rhs_ALLOC>
Matrix<ScalarType, device_name,  ALLOC>& Matrix<ScalarType, device_name,  ALLOC>::operator=(
    const Matrix<ScalarRhs, rhs_device_name, rhs_ALLOC>& rhs) {
  static_assert(sizeof(ScalarType) == sizeof(ScalarRhs),
                "sizeof ScalarType and ScalarRhs are not equal");
  resizeNoCopy(rhs.size_);
  util::memoryCopy(data_, leadingDimension(), rhs.data_, rhs.leadingDimension(), size_);
  return *this;
}

#endif

template <typename ScalarType, DeviceType device_name, class ALLOC>
bool Matrix<ScalarType, device_name, ALLOC>::operator==(const Matrix<ScalarType, device_name,  ALLOC>& other) const {
  if (device_name == GPU)
    return Matrix<ScalarType, CPU>(*this) == Matrix<ScalarType, CPU>(other);

  if (size() != other.size())
    return nrRows() * nrCols() == 0 and other.nrRows() * other.nrCols() == 0;

  for (int j = 0; j < nrCols(); ++j)
    for (int i = 0; i < nrRows(); ++i)
      if ((*this)(i, j) != other(i, j))
        return false;

  return true;
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
bool Matrix<ScalarType, device_name,  ALLOC>::operator!=(const Matrix<ScalarType, device_name,  ALLOC>& other) const {
  return not(*this == other);
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
void Matrix<ScalarType, device_name,  ALLOC>::resizeNoCopy(std::pair<int, int> new_size) {
  if (new_size.first > capacity_.first || new_size.second > capacity_.second) {
    size_ = new_size;
    capacity_ = capacityMultipleOfBlockSize(new_size);

    Allocator::deallocate(data_);
    data_ = Allocator::allocate(nrElements(capacity_));
  }
  else {
    size_ = new_size;
  }
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
void Matrix<ScalarType, device_name,  ALLOC>::clear() {
  Allocator::deallocate(data_);
  size_ = capacity_ = std::make_pair(0, 0);
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
void Matrix<ScalarType, device_name,  ALLOC>::swap(Matrix<ScalarType, device_name,  ALLOC>& rhs) {
  std::swap(size_, rhs.size_);
  std::swap(capacity_, rhs.capacity_);
  std::swap(data_, rhs.data_);
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
void Matrix<ScalarType, device_name,  ALLOC>::swapWithName(Matrix<ScalarType, device_name,  ALLOC>& rhs) {
  std::swap(name_, rhs.name_);
  swap(rhs);
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
template <DeviceType rhs_device_name, class rhs_ALLOC>
void Matrix<ScalarType, device_name,  ALLOC>::set(const Matrix<ScalarType, rhs_device_name, rhs_ALLOC>& rhs,
                                          int thread_id, int stream_id) {
  resize(rhs.size_);
  // This specialization is required since without unified memory CUDA doesn't known which memory locality the pointer has.
  if constexpr (device_name == DeviceType::GPU && rhs_device_name == DeviceType::CPU)
    util::memoryCopyH2D(data_, leadingDimension(), rhs.data_, rhs.leadingDimension(), size_,
                        thread_id, stream_id);
  else if constexpr (device_name == DeviceType::CPU && rhs_device_name == DeviceType::GPU)
    util::memoryCopyD2H(data_, leadingDimension(), rhs.data_, rhs.leadingDimension(), size_,
                        thread_id, stream_id);
  else if constexpr (device_name == DeviceType::CPU && rhs_device_name == DeviceType::CPU)
    util::memoryCopyCpu(data_, leadingDimension(), rhs.data_, rhs.leadingDimension(), size_,
                        thread_id, stream_id);
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
template <DeviceType rhs_device_name, class rhs_ALLOC>
void Matrix<ScalarType, device_name,  ALLOC>::set(const Matrix<ScalarType, rhs_device_name, rhs_ALLOC>& rhs,
                                          const util::GpuStream& stream [[maybe_unused]]) {
  resize(rhs.size_);
  if constexpr (device_name == DeviceType::GPU && rhs_device_name == DeviceType::CPU)
    util::memoryCopyH2D(data_, leadingDimension(), rhs.data_, rhs.leadingDimension(), size_);
  else if constexpr (device_name == DeviceType::CPU && rhs_device_name == DeviceType::GPU)
    util::memoryCopyD2H(data_, leadingDimension(), rhs.data_, rhs.leadingDimension(), size_);
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
template <DeviceType rhs_device_name, class rhs_ALLOC>
void Matrix<ScalarType, device_name,  ALLOC>::setAsync(const Matrix<ScalarType, rhs_device_name, rhs_ALLOC>& rhs,
                                               const util::GpuStream& stream) {
  resizeNoCopy(rhs.size_);
  util::memoryCopyAsync(data_, leadingDimension(), rhs.data_, rhs.leadingDimension(), size_, stream);
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
template <DeviceType rhs_device_name, class rhs_ALLOC>
void Matrix<ScalarType, device_name,  ALLOC>::setAsync(const Matrix<ScalarType, rhs_device_name, rhs_ALLOC>& rhs,
                                               const int thread_id, const int stream_id) {
  setAsync(rhs, util::getStream(thread_id, stream_id));
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
void Matrix<ScalarType, device_name,  ALLOC>::setToZero(const util::GpuStream& stream) {
  util::Memory<device_name>::setToZeroAsync(data_, leadingDimension() * nrCols(), stream);
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
void Matrix<ScalarType, device_name,  ALLOC>::print() const {
  if (device_name == GPU)
    return Matrix<ScalarType, CPU>(*this).print();

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

template <typename ScalarType, DeviceType device_name, class ALLOC>
std::string Matrix<ScalarType, device_name,  ALLOC>::toStr() const {
  if (device_name == GPU)
    return Matrix<ScalarType, CPU>(*this).toStr();

  std::stringstream ss;
  ss.precision(16);
  ss << std::scientific;

  ss << "\n";
  for (int i = 0; i < nrRows(); ++i) {
    for (int j = 0; j < nrCols(); ++j)
      ss << "\t" << operator()(i, j);
    ss << "\n";
  }

  return ss.str();
}
  
template <typename ScalarType, DeviceType device_name, class ALLOC>
void Matrix<ScalarType, device_name,  ALLOC>::printFingerprint() const {
  std::stringstream ss;

  ss << "\n";
  ss << "    name: " << name_ << "\n";
  ss << "    size: " << size_.first << ", " << size_.second << "\n";
  ss << "    capacity: " << capacity_.first << ", " << capacity_.second << "\n";
  ss << "    memory-size: " << nrElements(capacity_) * sizeof(ScalarType) * 1.e-6 << "(Mbytes)\n";

  std::cout << ss.str() << std::endl;
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
std::pair<int, int> Matrix<ScalarType, device_name,  ALLOC>::capacityMultipleOfBlockSize(
    std::pair<int, int> size) {
  assert(size.first >= 0);
  assert(size.second >= 0);

  auto get_new_size = [=](const int size) {
    return size <= 16 ? size : (size + block_size_ - 1) / block_size_ * block_size_;
  };

  size.first = get_new_size(size.first);
  size.second = get_new_size(size.second);

  return size;
}

template <typename ScalarType, DeviceType device_name, class ALLOC>
std::size_t Matrix<ScalarType, device_name,  ALLOC>::deviceFingerprint() const {
  if (device_name == GPU)
    return capacity_.first * capacity_.second * sizeof(ScalarType);
  else
    return 0;
}

/// Factory function for diangonal matrices, type is inferred from the type of Vector.
template <typename ScalarType, DeviceType device_name, class ALLOC>
auto makeDiagonalMatrix(Vector<ScalarType, device_name, ALLOC>& diag) {
  int dsize = diag.size();
  Matrix<ScalarType, device_name,  ALLOC> matrix("diag_matrix", dsize);
  for (int i = 0; i < dsize; ++i) {
    matrix(i, i) = diag[i];
  }
  return matrix;
}

/// Factory function for diangonal matrices, type is inferred from the type of Vector.
template <typename ScalarType, DeviceType device_name, class ALLOC>
auto makeDiagonalMatrixInv(Vector<ScalarType, device_name, ALLOC>& diag) {
  int dsize = diag.size();
  Matrix<ScalarType, device_name,  ALLOC> matrix("diag_matrix", dsize);
  // insure that if ScalarType is complex the 1 is as well.
  // then std::complex will give us a proper complex multiplicative inverse
  ScalarType the_one{};
  the_one += 1.0;
  for (int i = 0; i < dsize; ++i) {
    matrix(i, i) = the_one / diag[i];
  }
  return matrix;
}

}  // namespace linalg
}  // namespace dca

#endif  // DCA_LINALG_MATRIX_HPP
