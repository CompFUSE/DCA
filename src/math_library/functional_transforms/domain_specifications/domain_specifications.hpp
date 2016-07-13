// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_DOMAIN_SPECIFICATIONS_DOMAIN_SPECIFICATIONS_HPP
#define MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_DOMAIN_SPECIFICATIONS_DOMAIN_SPECIFICATIONS_HPP

#include <string>
#include <vector>

#include "math_library/functional_transforms/typedefs.hpp"

namespace math_algorithms {

template <typename scalartype, typename elementtype, DOMAIN_REPRESENTATIONS DMN_REP,
          BASIS_EXPANSIONS EXP, BOUNDARY_CONDITIONS BC, ELEMENT_SPACINGS SPACING>
class domain_specifications {
public:
  typedef scalartype scalar_type;
  typedef elementtype element_type;

  const static DOMAIN_REPRESENTATIONS DOMAIN_REPRESENTATION = DMN_REP;
  const static BASIS_EXPANSIONS BASIS_EXPANSION = EXP;

  const static BOUNDARY_CONDITIONS BOUNDARY_CONDITION = BC;
  const static ELEMENT_SPACINGS ELEMENT_SPACING = SPACING;

public:
  static std::string& to_str() {
    using math_algorithms::to_str;

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

}  // math_algorithms

#endif  // MATH_LIBRARY_FUNCTIONAL_TRANSFORMS_DOMAIN_SPECIFICATIONS_DOMAIN_SPECIFICATIONS_HPP
