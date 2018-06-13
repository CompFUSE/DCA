// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides methods to compute Gaussian probability functions.

#ifndef DCA_MATH_STATISTICS_GAUSSIAN_PROBABILITY_HPP
#define DCA_MATH_STATISTICS_GAUSSIAN_PROBABILITY_HPP

#include <array>
#include <cassert>
#include <cmath>

namespace dca {
namespace math {
namespace statistics {
namespace gauss {
namespace detail {
// dca::math::statistics::gauss::

inline double erf_inverse_help(const double t) {
  static const std::array<double, 3> c{2.515517, 0.802853, 0.010328};
  static const std::array<double, 3> d{1.432788, 0.189269, 0.001308};

  return t - ((c[2] * t + c[1]) * t + c[0]) / (((d[2] * t + d[1]) * t + d[0]) * t + 1.0);
}

inline double erf_inverse(const double p) {
  assert(p > 1.e-16 and p < 1 - 1.e-16);

  // See article below for explanation of this section.
  if (p < 0.5) {
    // F^-1(p) = - G^-1(p)
    return -erf_inverse_help(std::sqrt(-2.0 * std::log(p)));
  }
  else {
    // F^-1(p) = G^-1(1-p)
    return erf_inverse_help(std::sqrt(-2.0 * std::log(1 - p)));
  }
}

}  // detail
// dca::math::gauss::

// Finds x s.t. Q(x) = p where Q(x) = \frac{1}{\sqrt{2 \pi}} \int_x^\infty e^{-t^2/2} dt.
// Returns the result scaled for given mean 'mu' and standard deviation 'sigma'.
// The absolute value of the error of this rational approximation is smaller than 4.5 x 10^-4.
// See Abramowitz and Stegun, p. 933 eq. 26.2.23.
inline double argTailProbability(const double p, const double mu, const double sigma) {
  return sigma * detail::erf_inverse(p) + mu;
}

}  // gauss
}  // statistics
}  // math
}  // dca

#endif  // DCA_MATH_STATISTICS_GAUSSIAN_PROBABILITY_HPP
