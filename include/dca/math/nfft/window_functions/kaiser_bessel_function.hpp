// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// The Kaiser Bessel function is defined as,
// \f{eqnarray}{
// \phi(\tau) &=&
// \left\{
// \begin{array}{rl}
// \frac{1}{\pi\:\phi_0} \frac{\sinh(b\:\sqrt{m^2-(n\:\tau)^2})}{\sqrt{m^2-(n\:\tau)^2})} & \mbox{if
// } |\tau| < \frac{m}{n}
// \frac{1}{\pi\:\phi_0} \frac{\sin (b\:\sqrt{m^2-(n\:\tau)^2})}{\sqrt{m^2-(n\:\tau)^2})} & \mbox{if
// } |\tau| > \frac{m}{n}.
// \end{array} \right.
// \phi(\varpi) &=&
// \left\{
// \begin{array}{rl}
// \frac{1}{n\:\phi_0} I_0(m\:\sqrt{b^2-(2\:\pi\:k/n)^2}) & \mbox{for } k =
// -n\:(1-\frac{1}{2\sigma}), ..., n\:(1-\frac{1}{2\sigma})
// 0 & \mbox{ for all other k.}
// \end{array} \right.
// \f}
//
// Where  \f$I_0\f$ is the modified zero order Bessel function. Here, the symbols b and phi are,
//
// \f{eqnarray}{
// \phi_0 &=& \frac{1}{\pi\:m}\sinh\Big(\pi\:m\:\Big(\frac{2/\sigma-1}{\sigma}\Big)\Big)
// b &=& \pi\:\Big(2.-\frac{1}{\sigma}\Big)
// \f}

#ifndef DCA_MATH_NFFT_WINDOW_FUNCTIONS_KAISER_BESSEL_FUNCTION_HPP
#define DCA_MATH_NFFT_WINDOW_FUNCTIONS_KAISER_BESSEL_FUNCTION_HPP

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace dca {
namespace math {
namespace nfft {
// dca::math::nfft::

template <class Dnfft>
class kaiser_bessel_function {
public:
  static int n;
  static int m;
  static double sigma;

  static double phi_t(double x);
  static double d_phi_t(double x);
  static double phi_wn(int wn);

  static void test_besselI0();

private:
  static double renorm_factor() {
    return 1. / (M_PI * m) * std::sinh(M_PI * m * (2. * sigma - 1.) / sigma);
  }

  static double b_val() {
    return M_PI * (2. - 1. / sigma);
  }

  static double besselI0(double x);
};

template <class Dnfft>
int kaiser_bessel_function<Dnfft>::n = 1;
template <class Dnfft>
int kaiser_bessel_function<Dnfft>::m = 1;
template <class Dnfft>
double kaiser_bessel_function<Dnfft>::sigma = 1;

template <class Dnfft>
double kaiser_bessel_function<Dnfft>::phi_t(double x) {
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

template <class Dnfft>
double kaiser_bessel_function<Dnfft>::d_phi_t(double x) {
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

template <class Dnfft>
double kaiser_bessel_function<Dnfft>::phi_wn(int wn) {
  if (std::pow(wn, 2) < std::pow(n * (1. - 1. / (2. * sigma)), 2))
    return 1. / double(n) *
           besselI0(m * std::sqrt(std::pow(b_val(), 2) -
                                  std::pow((2. * M_PI * double(wn) / double(n)), 2))) /
           renorm_factor();
  else
    throw std::logic_error(__FUNCTION__);
}

template <class Dnfft>
double kaiser_bessel_function<Dnfft>::besselI0(double x) {
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

template <class Dnfft>
void kaiser_bessel_function<Dnfft>::test_besselI0() {
  for (int l = -m * b_val(); l <= m * b_val(); l++) {
    double x = l;
    std::cout << x << "\t" << besselI0(x) << "\n";
  }
}

}  // nfft
}  // math
}  // dca

#endif  // DCA_MATH_NFFT_WINDOW_FUNCTIONS_KAISER_BESSEL_FUNCTION_HPP
