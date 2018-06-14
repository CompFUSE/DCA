// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Union of piecewise linear and delta function.

#ifndef DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_UNION_PIECEWISE_LINEAR_AND_DELTA_FUNCTION_HPP
#define DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_UNION_PIECEWISE_LINEAR_AND_DELTA_FUNCTION_HPP

#include <complex>
#include <vector>

#include "dca/phys/dca_analysis/cpe_solver/basis_functions/delta_function.hpp"
#include "dca/phys/dca_analysis/cpe_solver/basis_functions/piecewise_linear_function.hpp"

namespace dca {
namespace phys {
namespace analysis {
// dca::phys::analysis::

class UnionPiecewiseLinearAndDeltaFunction {
public:
  using element_type = double;

  static int& get_size() {
    static int size = 0;
    return size;
  }

  static std::vector<double>& get_elements() {
    static std::vector<double> elements(0);
    return elements;
  }

  template <typename parameters_type>
  static void initialize(parameters_type& parameters);

  static double volume(int n);

  static std::complex<double> phi(int n, std::complex<double> z);
};

}  // analysis
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_UNION_PIECEWISE_LINEAR_AND_DELTA_FUNCTION_HPP
