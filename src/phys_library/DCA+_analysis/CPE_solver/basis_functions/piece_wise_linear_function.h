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

#ifndef PHYS_LIBRARY_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_PIECE_WISE_LINEAER_FUNCTION_H
#define PHYS_LIBRARY_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_PIECE_WISE_LINEAER_FUNCTION_H

#include <cassert>
#include <complex>
#include <vector>

#include "comp_library/function_library/domains/special_domains/dmn_0.h"
#include "phys_library/domains/time_and_frequency/frequency_domain_real_axis.h"

class piece_wise_linear_function {
public:
  using element_type = double;
  using w_REAL = dmn_0<frequency_domain_real_axis>;

public:
  static int& get_size();
  static std::vector<double>& get_elements();

  template <typename parameters_type>
  static void initialize(parameters_type& parameters);

  static double volume(int n);
  static std::complex<double> phi(int n, std::complex<double> z);
};

int& piece_wise_linear_function::get_size() {
  static int size = w_REAL::dmn_size();
  return size;
}

std::vector<double>& piece_wise_linear_function::get_elements() {
  static std::vector<double> elements = w_REAL::get_elements();
  return elements;
}

template <typename parameters_type>
void piece_wise_linear_function::initialize(parameters_type& /*parameters*/) {}

double piece_wise_linear_function::volume(int /*n*/) {
  // triangle with 2*delta_x/2.
  static double volume = 2. * (get_elements()[1] - get_elements()[0]) / 2.;
  return volume;
}

std::complex<double> piece_wise_linear_function::phi(int n, std::complex<double> z) {
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

#endif  // PHYS_LIBRARY_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_PIECE_WISE_LINEAER_FUNCTION_H
