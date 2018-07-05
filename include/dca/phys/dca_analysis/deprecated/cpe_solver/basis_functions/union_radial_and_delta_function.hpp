// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Union of radial and delta function.

#ifndef DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_UNION_RADIAL_AND_DELTA_FUNCTION_HPP
#define DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_UNION_RADIAL_AND_DELTA_FUNCTION_HPP

#include <complex>
#include <vector>

#include "dca/function/domains/dmn_0.hpp"
#include "dca/phys/dca_analysis/cpe_solver/basis_functions/delta_function.hpp"
#include "dca/phys/dca_analysis/cpe_solver/basis_functions/radial_function.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain_real_axis.hpp"

namespace dca {
namespace phys {
namespace analysis {
// dca::phys::analysis::

class UnionRadialAndDeltaFunction {
public:
  using element_type = double;
  using w_REAL = func::dmn_0<domains::frequency_domain_real_axis>;

  static int get_size() {
    return get_elements().size();
  }

  static std::vector<double>& get_elements() {
    static std::vector<double> elements(0);
    return elements;
  }

  template <typename parameters_type>
  static void initialize(parameters_type& parameters);

  static double volume(int n);

  static std::complex<double> phi(int n, std::complex<double> z);

private:
  static std::vector<double> poles;
};

template <typename parameters_type>
void UnionRadialAndDeltaFunction::initialize(parameters_type& parameters) {
  poles = parameters.get_poles();
  get_elements() = w_REAL::get_elements();
}

}  // analysis
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_UNION_RADIAL_AND_DELTA_FUNCTION_HPP
