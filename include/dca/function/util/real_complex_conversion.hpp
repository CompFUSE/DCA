// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// This file provides utility methods to convert functions from real to complex and vice versa.

#ifndef DCA_FUNCTION_UTIL_REAL_COMPLEX_CONVERSION_HPP
#define DCA_FUNCTION_UTIL_REAL_COMPLEX_CONVERSION_HPP

#include <complex>
#include <limits>
#include <stdexcept>

#include "dca/function/function.hpp"

namespace dca {
  enum class ImagCheck { IGNORE, WARN, FAIL };
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

namespace detail {
 /** implementation of function complex to real conversion function
  *  This allows good code reuse with no loss or low loss of performance
  *  in the ignore and warn cases.
  */
template <typename Scalartype, typename Dmn, ImagCheck IC>
auto real(const function<std::complex<Scalartype>, Dmn>& f) {
  function<Scalartype, Dmn> f_real;
  auto checkForImaginary = [](const std::complex<Scalartype> s) -> bool {
    return (std::abs(s.imag()) > (1000 * std::numeric_limits<Scalartype>::epsilon()));
  };
  auto writeCheckFail = [](const int i, const std::complex<Scalartype> s) -> std::string {
    std::ostringstream fail_message;
    fail_message << "Element " << i << " of function conversion complex -> real "
                 << " has a non negligible imaginary component of " << std::abs(s.imag()) << '\n';
    return fail_message.str();
  };
  bool have_not_warned = true;
  for (int i = 0; i < f_real.size(); ++i) {
    if constexpr (IC == ImagCheck::WARN) {
      if (have_not_warned) {
        if (checkForImaginary(f(i))) {
          std::cerr << "WARNING: " << writeCheckFail(i, f(i));
          have_not_warned = false;
        }
      }
    }
    else if constexpr (IC == ImagCheck::FAIL) {
      if (checkForImaginary(f(i)))
        throw(std::out_of_range(writeCheckFail(i, f(i))));
    }
    f_real(i) = f(i).real();
  }
  return f_real;
}
}  // namespace detail

// Returns a real valued function that is equal to the real part of f.
// If IC ==
//   ImagCheck::FAIL    -> throw exception if imaginary part of f is zero
//   ImagCheck::WARN    -> write one warning per conversion to std::err continue without further
//   checks ImagCheck::IGNORE  -> do not check for imaginary elements in f at all
template <typename Scalartype, typename Dmn>
auto real(const function<std::complex<Scalartype>, Dmn>& f, ImagCheck ic = ImagCheck::IGNORE) {
  switch (ic) {
    case ImagCheck::IGNORE:
      return detail::real<Scalartype, Dmn, ImagCheck::IGNORE>(f);
    case ImagCheck::WARN:
      return detail::real<Scalartype, Dmn, ImagCheck::WARN>(f);
    case ImagCheck::FAIL:
    default:
      return detail::real<Scalartype, Dmn, ImagCheck::FAIL>(f);
  }
}

}  // namespace util
}  // namespace func
}  // namespace dca

#endif  // DCA_FUNCTION_UTIL_REAL_COMPLEX_CONVERSION_HPP
