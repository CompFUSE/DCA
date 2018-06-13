// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements union_piecewise_linear_and_delta_function.hpp.

#include "dca/phys/dca_analysis/cpe_solver/basis_functions/union_piecewise_linear_and_delta_function.hpp"

#include <cassert>
#include <iostream>

namespace dca {
namespace phys {
namespace analysis {
// dca::phys::analysis::

template <typename parameters_type>
void UnionPiecewiseLinearAndDeltaFunction::initialize(parameters_type& parameters) {
  DeltaFunction::initialize(parameters);
  PiecewiseLinearFunction::initialize(parameters);

  get_size() = PiecewiseLinearFunction::get_size() + DeltaFunction::get_size();

  get_elements().insert(get_elements().end(), PiecewiseLinearFunction::get_elements().begin(),
                        PiecewiseLinearFunction::get_elements().end());
  get_elements().insert(get_elements().end(), DeltaFunction::get_elements().begin(),
                        DeltaFunction::get_elements().end());

  std::cout << __PRETTY_FUNCTION__ << "\n\t" << PiecewiseLinearFunction::get_size() << "\t"
            << DeltaFunction::get_size() << "\t" << get_size() << "\n";
}

double UnionPiecewiseLinearAndDeltaFunction::volume(int n) {
  assert(n >= 0 && n < get_size());
  return PiecewiseLinearFunction::volume(0);
}

std::complex<double> UnionPiecewiseLinearAndDeltaFunction::phi(int n, std::complex<double> z) {
  assert(n >= 0 && n < get_size());

  if (n < PiecewiseLinearFunction::get_size())
    return PiecewiseLinearFunction::phi(n, z);
  else
    return DeltaFunction::phi(n - PiecewiseLinearFunction::get_size(), z);
}

}  // analysis
}  // phys
}  // dca
