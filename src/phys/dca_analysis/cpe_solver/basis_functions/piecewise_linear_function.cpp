// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements piecewise_linear_function.hpp.

#include "dca/phys/dca_analysis/cpe_solver/basis_functions/piecewise_linear_function.hpp"

#include <cassert>

namespace dca {
namespace phys {
namespace analysis {
// dca::phys::analysis::

std::complex<double> PiecewiseLinearFunction::phi(int n, std::complex<double> z) {
  assert(n >= 0 && n < get_size());

  double delta_x = get_elements()[1] - get_elements()[0];

  std::complex<double> A_mn;

  std::complex<float> Z_fl(real(z), imag(z));

  std::complex<float> x0(w_REAL::get_elements()[n] - delta_x, 0.);
  std::complex<float> x1(w_REAL::get_elements()[n], 0.);
  std::complex<float> x2(w_REAL::get_elements()[n] + delta_x, 0.);

  std::complex<float> K_z =
      (-x0 + x1 + (x0 - Z_fl) * (std::log((x0 - Z_fl) / (x1 - Z_fl)))) / (x0 - x1) +
      (x1 - x2 - (x2 - Z_fl) * (std::log((-x1 + Z_fl) / (-x2 + Z_fl)))) / (x1 - x2);

  A_mn.real(real(K_z));
  A_mn.imag(imag(K_z));

  return A_mn;
}

}  // analysis
}  // phys
}  // dca
