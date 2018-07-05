// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Delta function.

#ifndef DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_DELTA_FUNCTION_HPP
#define DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_DELTA_FUNCTION_HPP

#include <complex>
#include <iostream>
#include <vector>

namespace dca {
namespace phys {
namespace analysis {
// dca::phys::analysis::

class DeltaFunction {
public:
  using element_type = double;

  static int& get_size() {
    static int size = 0;
    return size;
  }

  static std::vector<double>& get_elements() {
    static std::vector<double> elements(0, 0);
    return elements;
  }

  template <typename parameters_type>
  static void initialize(parameters_type& parameters);

  static std::complex<double> phi(int n, std::complex<double> z);

private:
  static double& epsilon() {
    static double epsilon = 1.e-3;
    return epsilon;
  }
};

template <typename parameters_type>
void DeltaFunction::initialize(parameters_type& parameters) {
  std::cout << __PRETTY_FUNCTION__ << std::endl;

  get_size() = parameters.get_poles().size();
  get_elements() = parameters.get_poles();

  epsilon() = 2. * parameters.get_EPSILON();
}

}  // analysis
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_DELTA_FUNCTION_HPP
