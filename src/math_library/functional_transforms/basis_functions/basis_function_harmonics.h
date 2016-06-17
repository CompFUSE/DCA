// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// Description

#ifndef MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_FUNCTIONS_BASIS_FUNCTION_HARMONICS_H
#define MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_FUNCTIONS_BASIS_FUNCTION_HARMONICS_H

#include <cmath>
#include <complex>
#include <vector>

#include "math_library/functional_transforms/basis_functions/basis_function_template.h"
#include "math_library/functional_transforms/typedefs.hpp"

namespace math_algorithms {
namespace functional_transforms {
// math_algorithms::functional_transforms::

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

public:
  static f_scalar_type execute(int i, int j) {
    return execute(lhs_dmn_type::get_elements()[i], rhs_dmn_type::get_elements()[j]);
  }

  static f_scalar_type execute(lhs_element_type& lh_elem, rhs_element_type& rh_elem) {
    const static f_scalar_type I(0, 1);

    f_scalar_type phase = dot_prod(lh_elem, rh_elem);

    return std::exp(I * phase);
  }

private:
  template <typename scalartype>
  inline static scalartype dot_prod(scalartype x, scalartype y) {
    return x * y;
  }

  template <typename scalartype>
  inline static scalartype dot_prod(std::vector<scalartype> x, std::vector<scalartype> y) {
    scalartype result = 0;
    for (size_t l = 0; l < x.size(); l++)
      result += x[l] * y[l];
    return result;
  }
};

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

public:
  static f_scalar_type execute(int i, int j) {
    return execute(lhs_dmn_type::get_elements()[i], rhs_dmn_type::get_elements()[j]);
  }

  static f_scalar_type execute(lhs_element_type& lh_elem, rhs_element_type& rh_elem) {
    const static f_scalar_type I(0, 1);

    f_scalar_type phase = dot_prod(lh_elem, rh_elem);

    // return std::exp(I*phase)/lhs_dmn_type::get_volume();
    return std::exp(-I * phase) / lhs_dmn_type::get_volume();
  }

private:
  template <typename scalartype>
  inline static scalartype dot_prod(scalartype x, scalartype y) {
    return x * y;
  }

  template <typename scalartype>
  inline static scalartype dot_prod(std::vector<scalartype> x, std::vector<scalartype> y) {
    scalartype result = 0;
    for (size_t l = 0; l < x.size(); l++)
      result += x[l] * y[l];
    return result;
  }
};
}  // functional_transforms
}  // math_algorithms

#endif  // MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_BASIS_FUNCTIONS_BASIS_FUNCTION_HARMONICS_H
