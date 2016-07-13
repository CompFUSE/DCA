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

#ifndef PHYS_LIBRARY_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_DELTA_FUNCTION_HPP
#define PHYS_LIBRARY_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_DELTA_FUNCTION_HPP

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

class delta_function {
public:
  using element_type = double;

public:
  static int& get_size();
  static std::vector<double>& get_elements();

  template <typename parameters_type>
  static void initialize(parameters_type& parameters);

  static double& epsilon();

  static std::complex<double> phi(int n, std::complex<double> z);
};

template <typename parameters_type>
void delta_function::initialize(parameters_type& parameters) {
  std::cout << __PRETTY_FUNCTION__ << std::endl;

  get_size() = parameters.get_poles().size();
  get_elements() = parameters.get_poles();

  epsilon() = 2. * parameters.get_EPSILON();
}

int& delta_function::get_size() {
  static int size = 0;
  return size;
}

std::vector<double>& delta_function::get_elements() {
  static std::vector<double> elements(0, 0);
  return elements;
}

double& delta_function::epsilon() {
  static double epsilon = 1.e-3;
  return epsilon;
}

std::complex<double> delta_function::phi(int n, std::complex<double> z) {
  assert(n >= 0 && n < get_size());

  double x = get_elements()[n];

  std::complex<double> A_mn;

  if (std::imag(z) > epsilon()) {
    A_mn.real(std::real(1. / (z - x)));
    A_mn.imag(std::imag(1. / (z - x)));
  }
  else {
    A_mn.real(std::fabs(std::real(z) - x) < 1.e-6 ? 0. : 1 / (std::real(z) - x));
    A_mn.imag(0.);
  }

  return A_mn;
}

#endif  // PHYS_LIBRARY_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_DELTA_FUNCTION_HPP
