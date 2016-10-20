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

#ifndef PHYS_LIBRARY_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_RADIAL_FUNCTION_HPP
#define PHYS_LIBRARY_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_RADIAL_FUNCTION_HPP

#include <cassert>
#include <cmath>
#include <complex>
#include <vector>

#include "dca/function/domains/dmn_0.hpp"
#include "dca/math/util/sgn.hpp"

#include "phys_library/domains/time_and_frequency/frequency_domain_real_axis.h"

class radial_function {
public:
  using element_type = double;
  using w_REAL = func::dmn_0<frequency_domain_real_axis>;

public:
  static int get_size();
  static std::vector<double>& get_elements();

  template <typename parameters_type>
  static void initialize(parameters_type& parameters);

  static double volume(int n);
  static std::complex<double> phi(int n, std::complex<double> z);
};

int radial_function::get_size() {
  return get_elements().size();
}

std::vector<double>& radial_function::get_elements() {
  static std::vector<double> elements(0, 0);
  return elements;
}

template <typename parameters_type>
void radial_function::initialize(parameters_type& parameters) {
  std::vector<double> elements = w_REAL::get_elements();

  for (size_t l = 0; l < elements.size(); l++) {
    elements[l] = math::util::sgn(elements[l]) *
                  std::pow(std::abs(elements[l]) / std::abs(elements[0]), 2) * std::abs(elements[0]);
  }

  get_elements() = elements;
}

double radial_function::volume(int n) {
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

std::complex<double> radial_function::phi(int n, std::complex<double> z) {
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

#endif  // PHYS_LIBRARY_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_RADIAL_FUNCTION_HPP
