// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides a struct that specifies domain representation, basis expansion, boundary
// condition, and element spacing.

#ifndef DCA_MATH_FUNCTION_TRANSFORM_DOMAIN_SPECIFICATIONS_HPP
#define DCA_MATH_FUNCTION_TRANSFORM_DOMAIN_SPECIFICATIONS_HPP

#include <string>
#include <utility>
#include <vector>

#include "dca/math/function_transform/basis_expansions.hpp"
#include "dca/math/function_transform/boundary_conditions.hpp"
#include "dca/math/function_transform/domain_representations.hpp"
#include "dca/math/function_transform/element_spacings.hpp"

namespace dca {
namespace math {
namespace transform {
// dca::math::transform::

template <typename scalartype, typename elementtype, DOMAIN_REPRESENTATIONS DMN_REP,
          BASIS_EXPANSIONS EXP, BOUNDARY_CONDITIONS BC, ELEMENT_SPACINGS SPACING>
struct domain_specifications {
  typedef scalartype scalar_type;
  typedef elementtype element_type;

  const static DOMAIN_REPRESENTATIONS DOMAIN_REPRESENTATION = DMN_REP;
  const static BASIS_EXPANSIONS BASIS_EXPANSION = EXP;

  const static BOUNDARY_CONDITIONS BOUNDARY_CONDITION = BC;
  const static ELEMENT_SPACINGS ELEMENT_SPACING = SPACING;

  static std::string& to_str() {
    using dca::math::transform::to_str;

    static std::string name = to_str(DOMAIN_REPRESENTATION) + " " + to_str(BASIS_EXPANSION) + " " +
                              to_str(BOUNDARY_CONDITION) + " " + to_str(ELEMENT_SPACING);

    return name;
  }
};

using interval_dmn_1D_type =
    domain_specifications<double, double, CONTINUOUS, HERMITE_CUBIC_SPLINE, INTERVAL, EQUIDISTANT>;
using interval_dmn_nD_type = domain_specifications<double, std::vector<double>, CONTINUOUS,
                                                   HERMITE_CUBIC_SPLINE, INTERVAL, EQUIDISTANT>;

using periodic_interval_dmn_1D_type =
    domain_specifications<double, double, CONTINUOUS, HERMITE_CUBIC_SPLINE, PERIODIC, EQUIDISTANT>;
using periodic_interval_dmn_nD_type =
    domain_specifications<double, std::vector<double>, CONTINUOUS, HERMITE_CUBIC_SPLINE, PERIODIC,
                          EQUIDISTANT>;

using harmonic_dmn_1D_type =
    domain_specifications<double, double, EXPANSION, HARMONICS, INTERVAL, EQUIDISTANT>;
using harmonic_dmn_nD_type =
    domain_specifications<double, std::vector<double>, EXPANSION, HARMONICS, INTERVAL, EQUIDISTANT>;

using legendre_dmn_1D_type =
    domain_specifications<double, double, EXPANSION, LEGENDRE_P, INTERVAL, EQUIDISTANT>;
using legendre_dmn_nD_type =
    domain_specifications<double, std::vector<double>, EXPANSION, LEGENDRE_P, INTERVAL, EQUIDISTANT>;

using Y_lm_dmn_1D_type =
    domain_specifications<double, std::pair<int, int>, EXPANSION, LEGENDRE_LM, INTERVAL, EQUIDISTANT>;
using Y_lm_dmn_nD_type = domain_specifications<double, std::vector<std::pair<int, int>>, EXPANSION,
                                               LEGENDRE_LM, INTERVAL, EQUIDISTANT>;

using discrete_interval_dmn_1D_type =
    domain_specifications<double, double, DISCRETE, KRONECKER_DELTA, INTERVAL, EQUIDISTANT>;
using discrete_interval_dmn_nD_type =
    domain_specifications<double, std::vector<double>, DISCRETE, KRONECKER_DELTA, INTERVAL, EQUIDISTANT>;

using discrete_periodic_dmn_1D_type =
    domain_specifications<double, double, DISCRETE, KRONECKER_DELTA, PERIODIC, EQUIDISTANT>;
using discrete_periodic_dmn_nD_type =
    domain_specifications<double, std::vector<double>, DISCRETE, KRONECKER_DELTA, PERIODIC, EQUIDISTANT>;

}  // transform
}  // math
}  // dca

#endif  // DCA_MATH_FUNCTION_TRANSFORM_DOMAIN_SPECIFICATIONS_HPP
