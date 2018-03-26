// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//
// This file provides the following helper methods based on the function class:
// - difference
// - complex
// - real

#ifndef DCA_FUNCTION_FUNCTION_UTILS_HPP
#define DCA_FUNCTION_FUNCTION_UTILS_HPP

#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>

#include "dca/function/function.hpp"

namespace dca {
namespace func {
namespace utils {
// dca::func::utils::

struct Difference {
  double l1, l2, l_inf;
};

// Returns the l1, l2, and l_inf norms of |f1 - f2| / |f1|.
template <typename Scalartype, class Dmn1, class Dmn2>
Difference difference(const function<Scalartype, Dmn1>& f1, const function<Scalartype, Dmn2>& f2) {
  if (f1.size() != f2.size())
    throw(std::logic_error("Sizes are different."));

  double l1 = 0.;
  double l2 = 0.;
  double linf = 0.;

  double l1_error = 0.;
  double l2_error = 0.;
  double linf_error = 0.;

  for (int ind = 0; ind < f2.size(); ind++) {
    const double f_abs = std::abs(f1(ind));
    l1 += f_abs;
    l2 += f_abs * f_abs;
    linf = std::max(linf, f_abs);

    const double err = std::abs(f1(ind) - f2(ind));
    l1_error += err;
    l2_error += err * err;
    linf_error = std::max(linf_error, err);
  }

  l1_error /= l1;
  l2_error = std::sqrt(l2_error / l2);
  linf_error /= linf;

  return Difference{l1_error, l2_error, linf_error};
}

// Returns a cast to complex of the function elements.
template <typename Scalartype, class Dmn>
function<std::complex<Scalartype>, Dmn> complex(const function<Scalartype, Dmn>& f) {
  function<std::complex<Scalartype>, Dmn> f_cmplx;

  for (int i = 0; i < f_cmplx.size(); ++i)
    f_cmplx(i) = std::complex<Scalartype>(f(i));

  return f_cmplx;
}

// Returns a real function whose elements are the real part of the elements of f.
// Throws a std::logic_error if check_imaginary is true, and any of the arguments has a non zero,
// up to rounding errors, imaginary part.
template <typename Scalartype, class Dmn>
function<Scalartype, Dmn> real(const function<std::complex<Scalartype>, Dmn>& f,
                               const bool check_imaginary = false) {
  function<Scalartype, Dmn> f_real;

  const Scalartype epsilon = 10 * std::numeric_limits<Scalartype>::epsilon();
  for (int i = 0; i < f_real.size(); ++i) {
    if (check_imaginary and std::abs(f(i).imag()) > epsilon)
      throw(std::logic_error("The function is not real."));
    f_real(i) = f(i).real();
  }

  return f_real;
}

}  // utils
}  // func
}  // dca

#endif  // DCA_FUNCTION_FUNCTION_UTILS_HPP
