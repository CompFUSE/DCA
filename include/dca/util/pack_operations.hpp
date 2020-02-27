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

// if_all<b1, b2, ...>::value is true only if all template arguments are true, otherwise false.
template <bool b1, bool... bs>
struct if_all {
  constexpr static bool value = b1 && if_all<bs...>::value;
};
template <bool b>
struct if_all<b> {
  constexpr static bool value = b;
};

// product(T1 a1, T2 a2, ...) returns the product of all its arguments. Equivalent to a1 * a2 * ...
template <typename T>
constexpr T product(T first) {
  return first;
}

template <typename T, class... Args>
constexpr T product(T first, Args... args) {
  return first * product<Args...>(args...);
}

// sum(T1 a1, T2 a2, ...) returns the sum of all its arguments. Equivalent to a1 + a2 + ...
template <typename T = void>
constexpr unsigned sum() {
  return 0;
}
template <typename T, class... Args>
constexpr T sum(T first, Args... args) {
  return first + sum<Args...>(args...);
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
