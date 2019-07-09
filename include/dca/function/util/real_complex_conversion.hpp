// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//
// This file provides utility methods to convert functions from real to complex and vice versa.

#ifndef DCA_FUNCTION_UTIL_REAL_COMPLEX_CONVERSION_HPP
#define DCA_FUNCTION_UTIL_REAL_COMPLEX_CONVERSION_HPP

#include <complex>
#include <limits>
#include <stdexcept>

#include "dca/function/function.hpp"

namespace dca {
namespace func {
namespace util {
// dca::func::util::

// Returns a complex valued function whose real part is equal to f.
template <typename Scalartype, typename Dmn>
auto complex(const function<Scalartype, Dmn>& f) {
  function<std::complex<Scalartype>, Dmn> f_complex;

  for (int i = 0; i < f_complex.size(); ++i)
    f_complex(i) = std::complex<Scalartype>(f(i), 0);

  return f_complex;
}

// Returns a real valued function that is equal to the real part of f.
// If check_imaginary = true, the method checks whether the imaginary part of f is zero and throws a
// std::logic_error if this is not the case.
template <typename Scalartype, typename Dmn>
auto real(const function<std::complex<Scalartype>, Dmn>& f, const bool check_imaginary = false) {
  function<Scalartype, Dmn> f_real;

  for (int i = 0; i < f_real.size(); ++i) {
    if (check_imaginary && std::abs(f(i).imag()) > 500 * std::numeric_limits<Scalartype>::epsilon())
      throw(std::logic_error("The function is not purely real."));

    f_real(i) = f(i).real();
  }

  return f_real;
}

}  // namespace util
}  // namespace func
}  // namespace dca

#endif  // DCA_FUNCTION_UTIL_REAL_COMPLEX_CONVERSION_HPP
