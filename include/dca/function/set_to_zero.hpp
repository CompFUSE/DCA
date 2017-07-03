// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
