// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides the following utilities for performing operations over a pack of arguments:
//
// - if_all
// - product

#ifndef DCA_UTIL_PACK_OPERATIONS_HPP
#define DCA_UTIL_PACK_OPERATIONS_HPP

#include "dca/util/type_list.hpp"

namespace dca {
namespace util {
// dca::util::

// returns is true only if all template arguments are true.
template <class... Args>
constexpr bool ifAll(Args... args) {
  return (args &&...);
}

// product(T1 a1, T2 a2, ...) returns the product of all its arguments. Equivalent to a1 * a2 * ...
template <class... Args>
constexpr auto product(Args... args) {
  return (args *...);
}

// sum(T1 a1, T2 a2, ...) returns the sum of all its arguments. Equivalent to a1 + a2 + ...
// sum() returns 0.
template <class... Args>
constexpr auto sum(Args... args) {
  return (0 + ... + args);
}

// size_sum.
// size_sum<T_1, ..., T_n> = \sum_{i=1}^n sizeof(T_i)
template <typename... Ts>
constexpr unsigned size_sum = sum(sizeof(Ts)...);

template <typename... Ts>
constexpr unsigned size_sum<mp_list<Ts...>> = sum(sizeof(Ts)...);

}  // namespace util
}  // namespace dca

#endif  // DCA_UTIL_PACK_OPERATIONS_HPP
