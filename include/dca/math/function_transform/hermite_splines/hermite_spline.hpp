// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// The files implements Hermite splines.

#ifndef DCA_MATH_FUNCTION_TRANSFORM_HERMITE_SPLINES_HERMITE_SPLINE_HPP
#define DCA_MATH_FUNCTION_TRANSFORM_HERMITE_SPLINES_HERMITE_SPLINE_HPP

#include <cmath>

namespace dca {
namespace math {
namespace transform {
// dca::math::transform::

class hermite_spline {
public:
  //  u( 0 < |s| < 1 ) = 1-s
  //  u( 1 < |s| )     = 0
  template <typename lh_scalar_type, typename rh_scalar_type>
  static lh_scalar_type linear(lh_scalar_type x, rh_scalar_type y, rh_scalar_type d) {
    rh_scalar_type delta = std::abs((y - x) / d);

    if (delta >= 0 - 1.e-6 and delta < 1)
      return 1 - delta;

    return 0;
  }

  // u( 0 < |s| < 1 ) = (a+2)*s^3-(a+3)*s^2+1
  // u( 1 < |s| < 2 ) = a*s^3-5as^2+8as-4a
  // u( 2 < |s| )     = 0
  //
  // f[s_, a_] := Piecewise[{{0 = (a + 2)*s^3 - (a + 3)*s^2 + 1, Abs[s] >= 2},
  //                         {a*s^3 - 5 a s^2 + 8 a s - 4 a    , 2 >= Abs[s] >= 1},
  // 	                     {(a + 2)*s^3 - (a + 3)*s^2 + 1    , 1 >= Abs[s] >= 0}}]
  template <typename lh_scalar_type, typename rh_scalar_type>
  static lh_scalar_type cubic(lh_scalar_type x, rh_scalar_type y, rh_scalar_type d, rh_scalar_type a) {
    rh_scalar_type delta = std::abs((y - x) / d);

    if (delta >= 0 - 1.e-6 and delta < 1)
      return (a + 2) * std::pow(delta, 3) - (a + 3) * std::pow(delta, 2) + 1;

    if (delta >= 1 and delta < 2)
      return a * std::pow(delta, 3) - 5 * a * std::pow(delta, 2) + 8 * a * std::pow(delta, 1) - 4 * a;

    return 0;
  }
};

}  // transform
}  // math
}  // dca

#endif  // DCA_MATH_FUNCTION_TRANSFORM_HERMITE_SPLINES_HERMITE_SPLINE_HPP
