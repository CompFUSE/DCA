// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Radial function.

#ifndef DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_RADIAL_FUNCTION_HPP
#define DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_RADIAL_FUNCTION_HPP

#include <cmath>
#include <complex>
#include <cstdlib>
#include <vector>

#include "dca/function/domains/dmn_0.hpp"
#include "dca/math/util/sgn.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain_real_axis.hpp"

namespace dca {
namespace phys {
namespace analysis {
// dca::phys::analysis::

class RadialFunction {
public:
  using element_type = double;
  using w_REAL = func::dmn_0<domains::frequency_domain_real_axis>;

  static int get_size() {
    return get_elements().size();
  }

  static std::vector<double>& get_elements() {
    static std::vector<double> elements(0, 0);
    return elements;
  }

  template <typename parameters_type>
  static void initialize(parameters_type& parameters);

  static double volume(int n);

  static std::complex<double> phi(int n, std::complex<double> z);
};

template <typename parameters_type>
void RadialFunction::initialize(parameters_type& parameters) {
  std::vector<double> elements = w_REAL::get_elements();

  for (std::size_t l = 0; l < elements.size(); l++) {
    elements[l] = math::util::sgn(elements[l]) *
                  std::pow(std::abs(elements[l]) / std::abs(elements[0]), 2) * std::abs(elements[0]);
  }

  get_elements() = elements;
}

}  // analysis
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_RADIAL_FUNCTION_HPP
