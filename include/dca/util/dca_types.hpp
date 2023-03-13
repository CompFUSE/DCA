// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter W.  Doak (doakpw@ornl.gov)
//

/** \file
 *  This file provides some "global" types for dca++
 *
 *  These type definitions should be such that they don't force a particular precision,
 *  complex/real, or model regime onto the code using them.
 */
#ifndef DCA_TYPE_HPP
#define DCA_TYPE_HPP

#include "dca/util/type_utils.hpp"

namespace dca {
namespace util {
// template <typename Scalar>
// using SignType = std::conditional_t<dca::util::IsComplex_t<Scalar>::value, Scalar, std::int8_t>;

template <typename T, typename = bool>
struct SignType_impl {};

template <typename T>
struct SignType_impl<T, IsComplex<T>> {
  using type = T;
};

template <typename T>
struct SignType_impl<T, IsReal<T>> {
  using type = std::int8_t;
};

template <typename T>
using SignType = typename SignType_impl<T>::type;

}  // namespace util
}  // namespace dca

#endif
