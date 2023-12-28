// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//

/** \file
 *  This file provides complex type utilities and definitions that must come
 *  after gpu_type_mapping
 */
#ifndef DCA_GPU_COMPLEX_TYPES_H
#define DCA_GPU_COMPLEX_TYPES_H

#include "dca/linalg/util/gpu_type_mapping.hpp"

#if defined(DCA_HAVE_CUDA)
#include "dca/linalg/util/complex_operators_cuda.cu.hpp"
#endif

namespace dca {
namespace util {

#ifdef DCA_HAVE_CUDA
template <typename T>
struct IsMagmaComplex_t : std::disjunction<std::is_same<magmaFloatComplex, T>,
                                           std::is_same<magmaDoubleComplex, T>, std::false_type> {};

template <typename T>
using IsMagmaComplex = std::enable_if_t<IsMagmaComplex_t<std::decay_t<T>>::value, bool>;

template <typename T>
struct TheOne<T, IsCUDAComplex<T>> {
  static constexpr T value{1.0, 0.0};
};

template <typename T>
struct TheZero<T, IsCUDAComplex<T>> {
  static constexpr T value{0.0, 0.0};
};

template <typename T>
std::enable_if_t<IsCUDAComplex_t<T>::value, void> makeOne(T& one) {
  one = T{1.0, 0.0};
}

template <typename T>
std::enable_if_t<IsCUDAComplex_t<T>::value, void> makeZero(T& zero) {
  zero = T{0.0, 0.0};
}

#endif

}  // namespace util
}  // namespace dca

#include "dca/util/type_help.hpp"
#endif
