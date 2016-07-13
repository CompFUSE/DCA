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

#ifndef PHYS_LIBRARY_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_UNION_PIECE_WISE_LINEAR_AND_DELTA_FUNCTION_H
#define PHYS_LIBRARY_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_UNION_PIECE_WISE_LINEAR_AND_DELTA_FUNCTIONH_

#include <cassert>
#include <complex>
#include <iostream>
#include <vector>

#include "phys_library/DCA+_analysis/CPE_solver/basis_functions/delta_function.hpp"
#include "phys_library/DCA+_analysis/CPE_solver/basis_functions/piece_wise_linear_function.h"

class union_piece_wise_linear_and_delta_function {
public:
  using element_type = double;

public:
  static int& get_size();
  static std::vector<double>& get_elements();

  template <typename parameters_type>
  static void initialize(parameters_type& parameters);

  static double volume(int n);

  static std::complex<double> phi(int n, std::complex<double> z);
};

int& union_piece_wise_linear_and_delta_function::get_size() {
  static int size = 0;
  return size;
}

std::vector<double>& union_piece_wise_linear_and_delta_function::get_elements() {
  static std::vector<double> elements(0);
  return elements;
}

template <typename parameters_type>
void union_piece_wise_linear_and_delta_function::initialize(parameters_type& parameters) {
  delta_function::initialize(parameters);
  piece_wise_linear_function::initialize(parameters);

  get_size() = piece_wise_linear_function::get_size() + delta_function::get_size();

  get_elements().insert(get_elements().end(), piece_wise_linear_function::get_elements().begin(),
                        piece_wise_linear_function::get_elements().end());
  get_elements().insert(get_elements().end(), delta_function::get_elements().begin(),
                        delta_function::get_elements().end());

  std::cout << __PRETTY_FUNCTION__ << "\n\t" << piece_wise_linear_function::get_size() << "\t"
       << delta_function::get_size() << "\t" << get_size() << "\n";
}

double union_piece_wise_linear_and_delta_function::volume(int n) {
  assert(n >= 0 && n < get_size());

  return piece_wise_linear_function::volume(0);
}

std::complex<double> union_piece_wise_linear_and_delta_function::phi(int n, std::complex<double> z) {
  assert(n >= 0 && n < get_size());

  if (n < piece_wise_linear_function::get_size())
    return piece_wise_linear_function::phi(n, z);
  else
    return delta_function::phi(n - piece_wise_linear_function::get_size(), z);
}

#endif  // PHYS_LIBRARY_DCA_ANALYSIS_CPE_SOLVER_BASIS_FUNCTIONS_UNION_PIECE_WISE_LINEAR_AND_DELTA_FUNCTION_H
