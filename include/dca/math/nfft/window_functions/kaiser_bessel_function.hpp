// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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

namespace dca {
namespace math {
namespace nfft {
// dca::math::nfft::

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

}  // nfft
}  // math
}  // dca

#endif  // DCA_MATH_NFFT_WINDOW_FUNCTIONS_KAISER_BESSEL_FUNCTION_HPP
