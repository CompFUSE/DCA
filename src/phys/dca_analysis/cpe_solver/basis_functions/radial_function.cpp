// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements radial_function.hpp.

#include "dca/phys/dca_analysis/cpe_solver/basis_functions/radial_function.hpp"

#include <cassert>
#include <cmath>

namespace dca {
namespace phys {
namespace analysis {
// dca::phys::analysis::

double RadialFunction::volume(int n) {
  assert(n >= 0 && n < w_REAL::dmn_size());

  double volume;

  if (n == 0)
    volume = 2. * (get_elements()[1] - get_elements()[0]) / 2.;
  else {
    if (n == get_size() - 1)
      volume = 2. * (get_elements()[n] - get_elements()[n - 1]) / 2.;
    else
      volume = (get_elements()[n + 1] - get_elements()[n - 1]) / 2.;
  }

  return volume;
}

std::complex<double> RadialFunction::phi(int n, std::complex<double> z) {
  assert(n >= 0 && n < get_size());

  std::complex<double> A_mn, x0, x1, x2;

  if (n == 0) {
    double delta_x = (get_elements()[1] - get_elements()[0]);

    x0 = get_elements()[0] - delta_x;
    x1 = get_elements()[0];
    x2 = get_elements()[0] + delta_x;
  }
  else {
    if (n == get_size() - 1) {
      double delta_x = (get_elements()[n] - get_elements()[n - 1]);

      x0 = get_elements()[n] - delta_x;
      x1 = get_elements()[n];
      x2 = get_elements()[n] + delta_x;
    }
    else {
      x0 = get_elements()[n - 1];
      x1 = get_elements()[n];
      x2 = get_elements()[n + 1];
    }
  }

  A_mn.real(std::real((-x0 + x1 + (x0 - z) * (std::log(x0 - z) - std::log(x1 - z))) / (x0 - x1) +
                      (x1 - x2 - (x2 - z) * (std::log(-x1 + z) - std::log(-x2 + z))) / (x1 - x2)));
  A_mn.imag(std::imag((-x0 + x1 + (x0 - z) * (std::log(x0 - z) - std::log(x1 - z))) / (x0 - x1) +
                      (x1 - x2 - (x2 - z) * (std::log(-x1 + z) - std::log(-x2 + z))) / (x1 - x2)));

  assert(A_mn == A_mn);  // no nan's !

  return A_mn;
}

}  // analysis
}  // phys
}  // dca
