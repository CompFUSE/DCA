// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class represents an AoS where each array has the same length but arbitrary type.

#ifndef DCA_LINALG_MULTI_VECTOR_HPP
#define DCA_LINALG_MULTI_VECTOR_HPP

#include "dca/linalg/vector.hpp"
#include "dca/linalg/util/cuda_stream.hpp"
#include "dca/util/type_list.hpp"
#include "dca/util/pack_operations.hpp"

namespace dca {
namespace linalg {
// dca::linalg::

template <DeviceType device, typename... Ts>
class MultiVector {
public:
  using Types = dca::util::Typelist<Ts...>;
  template <unsigned id>
  using Type = typename dca::util::TypeAt<id, Types>::type;

  // Initialize each sub-array with size n.
  MultiVector(std::size_t n = 0);

  // Resize the container so that each sub-array has size n, invalidating references and values.
  void resizeNoCopy(std::size_t n);

  // Copy the values of rhs asynchronously.
  template <DeviceType other_device>
  void setAsync(const MultiVector<other_device, Ts...>& rhs, const linalg::util::CudaStream& stream) {
    size_ = rhs.size_;
    data_.setAsync(rhs.data_, stream);
  }

  // Returns a pointer to the beginning of the id-th array
  // Preconditions: 0 <= id < length(Ts...).
  template <unsigned id>
  auto get() -> Type<id>*;
  template <unsigned id>
  auto get() const -> const Type<id>*;

  std::size_t size() const {
    return size_;
  }

  // Allows setAsync to access the data on another device.
  template <DeviceType other_device, typename... T2s>
  friend class MultiVector;

private:
  template <unsigned id>
  std::size_t offset() const;

  Vector<unsigned char, device> data_;
  std::size_t size_;
};

template <DeviceType device, typename... Ts>
MultiVector<device, Ts...>::MultiVector(std::size_t n) {
  resizeNoCopy(n);
}

template <DeviceType device, typename... Ts>
void MultiVector<device, Ts...>::resizeNoCopy(std::size_t n) {
  data_.resizeNoCopy(n * dca::util::size_sum<Ts...>);
  size_ = n;
}

template <DeviceType device, typename... Ts>
template <unsigned id>
auto MultiVector<device, Ts...>::get() -> Type<id>* {
  unsigned char* ptr = data_.ptr() + offset<id>();
  return reinterpret_cast<Type<id>*>(ptr);
}

template <DeviceType device, typename... Ts>
template <unsigned id>
auto MultiVector<device, Ts...>::get() const -> const Type<id>* {
  const unsigned char* ptr = data_.ptr() + offset<id>();
  return reinterpret_cast<const Type<id>*>(ptr);
}

template <DeviceType device, typename... Ts>
template <unsigned id>
std::size_t MultiVector<device, Ts...>::offset() const {
  static_assert(id < sizeof...(Ts), "Invalid sub-array id.");

  constexpr unsigned size_t_sum = dca::util::size_sum<dca::util::Sublist<id, Ts...>>;
  return size_ * size_t_sum;
}

}  // namespace linalg
}  // namespace dca

#endif  // DCA_LINALG_MULTI_VECTOR_HPP
