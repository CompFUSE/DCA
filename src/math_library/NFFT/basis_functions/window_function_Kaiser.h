//-*-C++-*-
// Author: Peter Staar

#ifndef MATH_LIBRARY_NFFT_BASIS_FUNCTIONS_WINDOW_FUNCTION_KAISER_H
#define MATH_LIBRARY_NFFT_BASIS_FUNCTIONS_WINDOW_FUNCTION_KAISER_H

#include <iostream>
#include <limits>
#include <stdexcept>
#include <cmath>

namespace math_algorithms {
namespace NFFT {
// math_algorithms::NFFT::

// The Kaiser Bessel function is defined as,
/*
  \f{eqnarray}{
  \phi(\tau) &=&
  \left\{
  \begin{array}{rl}
  \frac{1}{\pi\:\phi_0} \frac{\sinh(b\:\sqrt{m^2-(n\:\tau)^2})}{\sqrt{m^2-(n\:\tau)^2})} & \mbox{if
  } |\tau| < \frac{m}{n},\\
  \frac{1}{\pi\:\phi_0} \frac{\sin (b\:\sqrt{m^2-(n\:\tau)^2})}{\sqrt{m^2-(n\:\tau)^2})} & \mbox{if
  } |\tau| > \frac{m}{n}.
  \end{array} \right. \\
  \phi(\varpi) &=&
  \left\{
  \begin{array}{rl}
  \frac{1}{n\:\phi_0} I_0(m\:\sqrt{b^2-(2\:\pi\:k/n)^2}) & \mbox{for } k =
  -n\:(1-\frac{1}{2\sigma}), ..., n\:(1-\frac{1}{2\sigma})\\
  0 & \mbox{ for all other k.}
  \end{array} \right.
  \f}

  Where  \f$I_0\f$ is the modified zero order Bessel function. Here, the symbols b and phi are,

  \f{eqnarray}{
  \phi_0 &=& \frac{1}{\pi\:m}\sinh\Big(\pi\:m\:\Big(\frac{2/\sigma-1}{\sigma}\Big)\Big) \\
  b &=& \pi\:\Big(2.-\frac{1}{\sigma}\Big)
  \f}
 */

template <int DIMENSION>
struct kaiser_bessel_function {
public:
  static int n;
  static int m;

  static double sigma;

public:
  inline static double phi_t(double x);
  inline static double d_phi_t(double x);
  inline static double phi_wn(int wn);

  static void test_besselI0();

private:
  inline static double renorm_factor();
  inline static double b_val();

  static double besselI0(double x);
};

template <int DIMENSION>
int kaiser_bessel_function<DIMENSION>::n = 1;

template <int DIMENSION>
int kaiser_bessel_function<DIMENSION>::m = 1;

template <int DIMENSION>
double kaiser_bessel_function<DIMENSION>::sigma = 1;

template <int DIMENSION>
double kaiser_bessel_function<DIMENSION>::phi_t(double x) {
  double result(0);
  double val = square(m) - square(n * x);

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

template <int DIMENSION>
double kaiser_bessel_function<DIMENSION>::d_phi_t(double x) {
  double result(0);
  double val = square(m) - square(n * x);

  if (val > 1.e-12)
    result = 1. / M_PI * (n * n * x) * ((-b_val() * cosh(b_val() * sqrt(val))) / val +
                                        sinh(b_val() * sqrt(val)) / (std::pow(val, 3. / 2.)));

  if (val < -1.e-12)
    result = 1. / M_PI * (n * n * x) * ((b_val() * cos(b_val() * sqrt(-val))) / (-val) -
                                        sin(b_val() * sqrt(-val)) / (std::pow(-val, 3. / 2.)));

  if (val >= -1.e-12 && val < 1.e-12)
    return 0;

  if (!(result == result))
    throw std::logic_error(__FUNCTION__);

  if (std::numeric_limits<double>::has_infinity && result == std::numeric_limits<double>::infinity())
    throw std::logic_error(__FUNCTION__);

  return result / renorm_factor();
}

template <int DIMENSION>
double kaiser_bessel_function<DIMENSION>::phi_wn(int wn) {
  if (square(wn) < square(n * (1. - 1. / (2. * sigma))))
    return 1. / double(n) *
           besselI0(m * sqrt(square(b_val()) - square((2. * M_PI * double(wn) / double(n))))) /
           renorm_factor();
  else
    throw std::logic_error(__FUNCTION__);
}

template <int DIMENSION>
double kaiser_bessel_function<DIMENSION>::b_val() {
  return M_PI * (2. - 1. / sigma);
}

template <int DIMENSION>
double kaiser_bessel_function<DIMENSION>::renorm_factor() {
  return 1. / (M_PI * m) * sinh(M_PI * m * (2. * sigma - 1.) / sigma);
}

template <int DIMENSION>
double kaiser_bessel_function<DIMENSION>::besselI0(double x) {
  double result = 0.;

  double sq_x_div_two = square(x / 2.);

  double gamma = 1;
  double y = 1.;

  double term;

  for (int l = 0; l < 1.e16; l++) {
    term = y / gamma;

    if (std::fabs(term) < 1.e-16)
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

template <int DIMENSION>
void kaiser_bessel_function<DIMENSION>::test_besselI0() {
  for (int l = -m * b_val(); l <= m * b_val(); l++) {
    double x = l;
    std::cout << x << "\t" << besselI0(x) << "\n";
  }
}

}  // NFFT
}  // math_algorithm

#endif  // MATH_LIBRARY_NFFT_BASIS_FUNCTIONS_WINDOW_FUNCTION_KAISER_H
