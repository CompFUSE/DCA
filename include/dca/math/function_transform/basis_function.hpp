// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides different types of basis functions.

#ifndef DCA_MATH_FUNCTION_TRANSFORM_BASIS_FUNCTION_HPP
#define DCA_MATH_FUNCTION_TRANSFORM_BASIS_FUNCTION_HPP

#include <cassert>
#include <cmath>
#include <complex>
#include <vector>

#include "dca/math/function_transform/basis_expansions.hpp"
#include "dca/math/function_transform/hermite_splines/hermite_cubic_spline.hpp"
#include "dca/math/function_transform/hermite_splines/hermite_spline.hpp"
#include "dca/math/util/vector_operations.hpp"

namespace dca {
namespace math {
namespace transform {
// dca::math::transform::

//
// Empty class template.
//
template <typename lhs_dmn_type, BASIS_EXPANSIONS BS_LHS, typename rhs_dmn_type, BASIS_EXPANSIONS BS_RHS>
class basis_function {};

//
// Template specialization for Legendre_P.
//
template <typename lh_dmn_type, BASIS_EXPANSIONS BS_LHS, typename rh_dmn_type>
class basis_function<lh_dmn_type, BS_LHS, rh_dmn_type, LEGENDRE_P> {
public:
  typedef typename lh_dmn_type::dmn_specifications_type lh_spec_dmn_type;
  typedef typename rh_dmn_type::dmn_specifications_type rh_spec_dmn_type;

  typedef typename lh_spec_dmn_type::scalar_type lh_scalar_type;
  typedef typename rh_spec_dmn_type::scalar_type rh_scalar_type;

  typedef typename lh_spec_dmn_type::element_type lh_element_type;
  typedef typename rh_spec_dmn_type::element_type rh_element_type;

  typedef lh_scalar_type f_scalar_type;

  static f_scalar_type execute(int i, int j);
  static f_scalar_type execute(lh_element_type& lh_elem, rh_element_type& rh_elem);

private:
  static f_scalar_type compute(lh_scalar_type t, int l);
};

//
// Template specialization for Hermite cubic splines.
//
template <typename lhs_dmn_type, BASIS_EXPANSIONS BS_LHS, typename rhs_dmn_type>
class basis_function<lhs_dmn_type, BS_LHS, rhs_dmn_type, HERMITE_CUBIC_SPLINE> {
public:
  typedef typename lhs_dmn_type::dmn_specifications_type lhs_spec_dmn_type;
  typedef typename rhs_dmn_type::dmn_specifications_type rhs_spec_dmn_type;

  typedef typename lhs_spec_dmn_type::scalar_type lhs_scalar_type;
  typedef typename rhs_spec_dmn_type::scalar_type rhs_scalar_type;

  typedef typename lhs_spec_dmn_type::element_type lhs_element_type;
  typedef typename rhs_spec_dmn_type::element_type rhs_element_type;

  const static ELEMENT_SPACINGS ELEMENT_SPACING = rhs_spec_dmn_type::ELEMENT_SPACING;
  const static BOUNDARY_CONDITIONS BOUNDARY_CONDITION = rhs_spec_dmn_type::BOUNDARY_CONDITION;

  typedef rhs_scalar_type f_scalar_type;

  static rhs_scalar_type execute(int i, int j);
  static rhs_scalar_type execute(lhs_element_type& x, rhs_element_type& y);
};

//
// Template specialization for Kronecker delta and harmonics.
//
template <typename lhs_dmn_type, typename rhs_dmn_type>
class basis_function<lhs_dmn_type, KRONECKER_DELTA, rhs_dmn_type, HARMONICS> {
public:
  typedef typename lhs_dmn_type::dmn_specifications_type lhs_spec_dmn_type;
  typedef typename rhs_dmn_type::dmn_specifications_type rhs_spec_dmn_type;

  typedef typename lhs_spec_dmn_type::scalar_type lhs_scalar_type;
  typedef typename rhs_spec_dmn_type::scalar_type rhs_scalar_type;

  typedef typename lhs_spec_dmn_type::element_type lhs_element_type;
  typedef typename rhs_spec_dmn_type::element_type rhs_element_type;

  typedef std::complex<lhs_scalar_type> f_scalar_type;

  static f_scalar_type execute(int i, int j);
  static f_scalar_type execute(lhs_element_type& lh_elem, rhs_element_type& rh_elem);
};

//
// Template specialization for Hermite cubic splines and harmonics.
//
template <typename lhs_dmn_type, typename rhs_dmn_type>
class basis_function<lhs_dmn_type, HERMITE_CUBIC_SPLINE, rhs_dmn_type, HARMONICS> {
public:
  typedef typename lhs_dmn_type::dmn_specifications_type lhs_spec_dmn_type;
  typedef typename rhs_dmn_type::dmn_specifications_type rhs_spec_dmn_type;

  typedef typename lhs_spec_dmn_type::scalar_type lhs_scalar_type;
  typedef typename rhs_spec_dmn_type::scalar_type rhs_scalar_type;

  typedef typename lhs_spec_dmn_type::element_type lhs_element_type;
  typedef typename rhs_spec_dmn_type::element_type rhs_element_type;

  typedef std::complex<lhs_scalar_type> f_scalar_type;

  static f_scalar_type execute(int i, int j);

  static f_scalar_type execute(lhs_element_type& lh_elem, rhs_element_type& rh_elem);
};

//
// Definitions of specialization for Legendre_P.
//
template <typename lh_dmn_type, BASIS_EXPANSIONS BS_LHS, typename rh_dmn_type>
typename basis_function<lh_dmn_type, BS_LHS, rh_dmn_type, LEGENDRE_P>::f_scalar_type basis_function<
    lh_dmn_type, BS_LHS, rh_dmn_type, LEGENDRE_P>::execute(int i, int j) {
  return execute(lh_dmn_type::get_elements()[i], rh_dmn_type::get_elements()[j]);
}

template <typename lh_dmn_type, BASIS_EXPANSIONS BS_LHS, typename rh_dmn_type>
typename basis_function<lh_dmn_type, BS_LHS, rh_dmn_type, LEGENDRE_P>::f_scalar_type basis_function<
    lh_dmn_type, BS_LHS, rh_dmn_type, LEGENDRE_P>::execute(lh_element_type& lh_elem,
                                                           rh_element_type& rh_elem) {
  lh_scalar_type delta = lh_dmn_type::get_volume();
  lh_scalar_type min = lh_dmn_type::get_min()[0];

  return compute(2 * (lh_elem - min) / delta - 1, rh_elem);
}

template <typename lh_dmn_type, BASIS_EXPANSIONS BS_LHS, typename rh_dmn_type>
typename basis_function<lh_dmn_type, BS_LHS, rh_dmn_type, LEGENDRE_P>::f_scalar_type basis_function<
    lh_dmn_type, BS_LHS, rh_dmn_type, LEGENDRE_P>::compute(lh_scalar_type t, int l) {
  assert(t > -1. - 1.e-6 and t < 1. + 1.e-6);

  f_scalar_type f_t = 0;
  lh_scalar_type l_v = l;

  switch (l) {
    case 0:
      f_t = 1;
      break;

    case 1:
      f_t = t;
      break;

    default:
      f_t += ((2 * l_v - 1) * t * compute(t, l - 1) - (l_v - 1) * compute(t, l - 2)) / l_v;
  }

  return f_t;
}

//
// Definitions of specialization for Hermite cubic splines.
//
template <typename lhs_dmn_type, BASIS_EXPANSIONS BS_LHS, typename rhs_dmn_type>
typename basis_function<lhs_dmn_type, BS_LHS, rhs_dmn_type, HERMITE_CUBIC_SPLINE>::rhs_scalar_type basis_function<
    lhs_dmn_type, BS_LHS, rhs_dmn_type, HERMITE_CUBIC_SPLINE>::execute(int i, int j) {
  return hermite_cubic_spline<lhs_dmn_type, rhs_dmn_type, ELEMENT_SPACING, BOUNDARY_CONDITION,
                              rhs_dmn_type::DIMENSION>::execute(i, j);
}

template <typename lhs_dmn_type, BASIS_EXPANSIONS BS_LHS, typename rhs_dmn_type>
typename basis_function<lhs_dmn_type, BS_LHS, rhs_dmn_type, HERMITE_CUBIC_SPLINE>::rhs_scalar_type basis_function<
    lhs_dmn_type, BS_LHS, rhs_dmn_type, HERMITE_CUBIC_SPLINE>::execute(lhs_element_type& x,
                                                                       rhs_element_type& y) {
  const static lhs_scalar_type a = -0.5;

  lhs_scalar_type d = rhs_dmn_type::get_elements()[1] - rhs_dmn_type::get_elements()[0];

  return hermite_spline::cubic(x, y, d, a);
}

//
// Definitions of specialization for Kronecker delta and harmonics.
//
template <typename lhs_dmn_type, typename rhs_dmn_type>
typename basis_function<lhs_dmn_type, KRONECKER_DELTA, rhs_dmn_type, HARMONICS>::f_scalar_type basis_function<
    lhs_dmn_type, KRONECKER_DELTA, rhs_dmn_type, HARMONICS>::execute(int i, int j) {
  return execute(lhs_dmn_type::get_elements()[i], rhs_dmn_type::get_elements()[j]);
}

template <typename lhs_dmn_type, typename rhs_dmn_type>
typename basis_function<lhs_dmn_type, KRONECKER_DELTA, rhs_dmn_type, HARMONICS>::f_scalar_type basis_function<
    lhs_dmn_type, KRONECKER_DELTA, rhs_dmn_type, HARMONICS>::execute(lhs_element_type& lh_elem,
                                                                     rhs_element_type& rh_elem) {
  const static f_scalar_type I(0, 1);

  f_scalar_type phase = util::innerProduct(lh_elem, rh_elem);

  return std::exp(I * phase);
}

//
// Definitions of specialization for Hermite cubic splines and harmonics.
//
template <typename lhs_dmn_type, typename rhs_dmn_type>
typename basis_function<lhs_dmn_type, HERMITE_CUBIC_SPLINE, rhs_dmn_type, HARMONICS>::f_scalar_type basis_function<
    lhs_dmn_type, HERMITE_CUBIC_SPLINE, rhs_dmn_type, HARMONICS>::execute(int i, int j) {
  return execute(lhs_dmn_type::get_elements()[i], rhs_dmn_type::get_elements()[j]);
}

template <typename lhs_dmn_type, typename rhs_dmn_type>
typename basis_function<lhs_dmn_type, HERMITE_CUBIC_SPLINE, rhs_dmn_type, HARMONICS>::f_scalar_type basis_function<
    lhs_dmn_type, HERMITE_CUBIC_SPLINE, rhs_dmn_type, HARMONICS>::execute(lhs_element_type& lh_elem,
                                                                          rhs_element_type& rh_elem) {
  const static f_scalar_type I(0, 1);

  f_scalar_type phase = util::innerProduct(lh_elem, rh_elem);

  // return std::exp(I*phase)/lhs_dmn_type::get_volume();
  return std::exp(-I * phase) / lhs_dmn_type::get_volume();
}

}  // transform
}  // math
}  // dca

#endif  // DCA_MATH_FUNCTION_TRANSFORM_BASIS_FUNCTION_HPP
