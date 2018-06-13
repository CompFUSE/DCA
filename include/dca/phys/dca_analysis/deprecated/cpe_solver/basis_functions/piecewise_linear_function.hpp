// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Piecewice linear function.

#ifndef DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_PIECEWISE_LINEAR_FUNCTION_HPP
#define DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_PIECEWISE_LINEAR_FUNCTION_HPP

#include <complex>
#include <vector>

#include "dca/function/domains/dmn_0.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain_real_axis.hpp"

namespace dca {
namespace phys {
namespace analysis {
// dca::phys::analysis::

class PiecewiseLinearFunction {
public:
  using element_type = double;
  using w_REAL = func::dmn_0<domains::frequency_domain_real_axis>;

  static int& get_size() {
    static int size = w_REAL::dmn_size();
    return size;
  }

  static std::vector<double>& get_elements() {
    static std::vector<double> elements = w_REAL::get_elements();
    return elements;
  }

  static double volume(int n) {
    // triangle with 2*delta_x/2.
    static double volume = 2. * (get_elements()[1] - get_elements()[0]) / 2.;
    return volume;
  }

  template <typename parameters_type>
  static void initialize(parameters_type& /*parameters*/) {}

  static std::complex<double> phi(int n, std::complex<double> z);
};

}  // analysis
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_PIECEWISE_LINEAR_FUNCTION_HPP
