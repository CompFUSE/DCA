// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements delta_function.hpp.

#include "dca/phys/dca_analysis/cpe_solver/basis_functions/delta_function.hpp"

#include <cassert>
#include <cmath>

namespace dca {
namespace phys {
namespace analysis {
// dca::phys::analysis::

std::complex<double> DeltaFunction::phi(int n, std::complex<double> z) {
  assert(n >= 0 && n < get_size());

  double x = get_elements()[n];

  std::complex<double> A_mn;

  if (std::imag(z) > epsilon()) {
    A_mn.real(std::real(1. / (z - x)));
    A_mn.imag(std::imag(1. / (z - x)));
  }
  else {
    A_mn.real(std::abs(std::real(z) - x) < 1.e-6 ? 0. : 1 / (std::real(z) - x));
    A_mn.imag(0.);
  }

  return A_mn;
}

}  // analysis
}  // phys
}  // dca
