// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements kaiser_bessel_function.hpp

#include "dca/math/nfft/window_functions/kaiser_bessel_function.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace dca {
namespace math {
namespace nfft {
// dca::math::nfft

int kaiser_bessel_function::n = 1;
int kaiser_bessel_function::m = 1;
double kaiser_bessel_function::sigma = 1;

double kaiser_bessel_function::phi_t(double x) {
  double result(0);
  double val = m * m - n * n * x * x;

  if (val > 1.e-12)
    result = 1. / M_PI * std::sinh(b_val() * std::sqrt(val)) / std::sqrt(val);

  if (val < -1.e-12)
    result = 1. / M_PI * std::sin(b_val() * std::sqrt(-val)) / std::sqrt(-val);

  if (val >= -1.e-12 && val < 1.e-12)
    result = b_val() / M_PI;

  if (!(result == result))
    throw std::logic_error(__FUNCTION__);

  if (std::numeric_limits<double>::has_infinity && result == std::numeric_limits<double>::infinity())
    throw std::logic_error(__FUNCTION__);

  return result / renorm_factor();
}

double kaiser_bessel_function::d_phi_t(double x) {
  double result(0);
  double val = m * m - n * n * x * x;

  if (val > 1.e-12)
    result =
        1. / M_PI * (n * n * x) * ((-b_val() * std::cosh(b_val() * std::sqrt(val))) / val +
                                   std::sinh(b_val() * std::sqrt(val)) / (std::pow(val, 3. / 2.)));

  if (val < -1.e-12)
    result =
        1. / M_PI * (n * n * x) * ((b_val() * std::cos(b_val() * std::sqrt(-val))) / (-val) -
                                   std::sin(b_val() * std::sqrt(-val)) / (std::pow(-val, 3. / 2.)));

  if (val >= -1.e-12 && val < 1.e-12)
    return 0;

  if (!(result == result))
    throw std::logic_error(__FUNCTION__);

  if (std::numeric_limits<double>::has_infinity && result == std::numeric_limits<double>::infinity())
    throw std::logic_error(__FUNCTION__);

  return result / renorm_factor();
}

double kaiser_bessel_function::phi_wn(int wn) {
  if (std::pow(wn, 2) < std::pow(n * (1. - 1. / (2. * sigma)), 2))
    return 1. / double(n) *
           besselI0(m * std::sqrt(std::pow(b_val(), 2) -
                                  std::pow((2. * M_PI * double(wn) / double(n)), 2))) /
           renorm_factor();
  else
    throw std::logic_error(__FUNCTION__);
}

double kaiser_bessel_function::besselI0(double x) {
  double result = 0.;

  double sq_x_div_two = x * x / (2. * 2.);

  double gamma = 1;
  double y = 1.;

  double term;

  for (int l = 0; l < 1.e16; l++) {
    term = y / gamma;

    if (std::abs(term) < 1.e-16)
      break;

    result += term;

    gamma *= ((l + 1) * (l + 1));
    y *= sq_x_div_two;
  }

  if (!(result == result))
    throw std::logic_error(__FUNCTION__);

  if (std::numeric_limits<double>::has_infinity && result == std::numeric_limits<double>::infinity())
    throw std::logic_error(__FUNCTION__);

  return result;
}

void kaiser_bessel_function::test_besselI0() {
  for (int l = -m * b_val(); l <= m * b_val(); l++) {
    double x = l;
    std::cout << x << "\t" << besselI0(x) << "\n";
  }
}

}  // nfft
}  // math
}  // dca
