// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Bart Ydens
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class assures that all functions with a scalar or complex type are initialized to zero.
// Functions with a different type are not initialized.
//
// INTERNAL: Require operator=(const scalartype&) for function's template parameter scalartype and
//           remove this file?

#ifndef DCA_FUNCTION_SET_TO_ZERO_HPP
#define DCA_FUNCTION_SET_TO_ZERO_HPP

#include <complex>
#include <type_traits>

namespace dca {
namespace func {
// dca::func::

template <class T>
inline void setToZero(std::complex<T>& whatever) {
  whatever = T(0);
}

template <class T>
inline std::enable_if_t<std::is_scalar<T>::value, void> setToZero(T& whatever) {
  whatever = T(0);
}

template <class T>
inline std::enable_if_t<not std::is_scalar<T>::value, void> setToZero(T& /*whatever*/) {}

}  // func
}  // dca

#endif  // DCA_FUNCTION_SET_TO_ZERO_HPP
