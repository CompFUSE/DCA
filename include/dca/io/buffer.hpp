// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides specific tests for the HDF5 reader.

#ifndef DCA_IO_BUFFER_HPP
#define DCA_IO_BUFFER_HPP

#include <algorithm>
#include <memory>
#include <vector>

namespace dca {
namespace io {
// dca::io::

class Buffer : public std::vector<unsigned char> {
public:
  using Container = std::vector<unsigned char>;

  // Copies obj in the buffer as a plain sequence of bytes.
  template <class T, typename = std::enable_if_t<std::is_trivially_copyable<T>::value>>
  Buffer& operator<<(const T& obj);

  template <class T>
  Buffer& operator<<(const std::unique_ptr<T>& ptr);
  template <class T>
  Buffer& operator<<(const std::vector<T>& vec);
  template <class T1, class T2>
  Buffer& operator<<(const std::pair<T1, T2>& p);

  // Reads obj from the buffer as a plain sequence of bytes.
  template <class T, typename = std::enable_if_t<std::is_trivially_copyable<T>::value>>
  Buffer& operator>>(T& obj);

  template <class T>
  Buffer& operator>>(std::unique_ptr<T>& ptr);
  template <class T>
  Buffer& operator>>(std::vector<T>& vec);
  template <class T1, class T2>
  Buffer& operator>>(std::pair<T1, T2>& p);

  // Retrieves the position of the next read.
  std::size_t tellg() const {
    return read_idx_;
  }
  // Sets the position of the next read.
  void setg(const std::size_t idx) {
    read_idx_ = idx;
  }

private:
  // Copies n objects of type T in the buffer as a plain sequence of bytes.
  template <class T>
  void write(const T* ptr, std::size_t n = 1);
  // Read n objects o f type T from the buffer as a plain sequence of bytes.
  template <class T>
  void read(T* ptr, std::size_t n = 1);

  std::size_t read_idx_ = 0;
};

template <class T>
void Buffer::write(const T* ptr, const std::size_t n) {
  const auto old_size = size();
  resize(old_size + sizeof(T) * n);
  std::copy_n(reinterpret_cast<const unsigned char*>(ptr), sizeof(T) * n, data() + old_size);
}

template <class T, typename>
Buffer& Buffer::operator<<(const T& obj) {
  write(&obj);
  return *this;
}

template <class T>
Buffer& Buffer::operator<<(const std::unique_ptr<T>& ptr) {
  const bool is_set(ptr);
  write(&is_set);
  if (is_set)
    *this << *ptr;
  return *this;
}

template <class T>
Buffer& Buffer::operator<<(const std::vector<T>& vec) {
  *this << vec.size();
  if (std::is_trivially_copyable<T>::value)
    write(vec.data(), vec.size());
  else
    for (const auto& elem : vec)
      *this << elem;

  return *this;
}

template <class T1, class T2>
Buffer& Buffer::operator<<(const std::pair<T1, T2>& p) {
  return (*this) << p.first << p.second;
}

template <class T>
void Buffer::read(T* ptr, const std::size_t n) {
  const std::size_t read_size = sizeof(T) * n;
  if (read_idx_ + read_size > size())
    throw(std::out_of_range("Buffer has no more data to read."));

  std::copy_n(data() + read_idx_, read_size, reinterpret_cast<unsigned char*>(ptr));
  read_idx_ += read_size;
}

template <class T, typename>
Buffer& Buffer::operator>>(T& obj) {
  read(&obj);
  return *this;
}

template <class T>
Buffer& Buffer::operator>>(std::unique_ptr<T>& ptr) {
  bool is_set;
  read(&is_set);
  if (is_set) {
    ptr.reset(new T);
    *this >> *ptr;
  }
  return *this;
}

template <class T>
Buffer& Buffer::operator>>(std::vector<T>& vec) {
  std::size_t vec_size;
  *this >> vec_size;
  vec.resize(vec_size);

  if (std::is_trivially_copyable<T>::value)
    read(vec.data(), vec_size);
  else
    for (auto& elem : vec)
      *this >> elem;

  return *this;
}

template <class T1, class T2>
Buffer& Buffer::operator>>(std::pair<T1, T2>& p) {
  return (*this) >> p.first >> p.second;
}

}  // io
}  // dca

#endif  // DCA_IO_BUFFER_HPP
