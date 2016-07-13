// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// Description

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ADVANCED_ED_STRUCTURES_PSI_SYMMETRY_OPERATOR_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ADVANCED_ED_STRUCTURES_PSI_SYMMETRY_OPERATOR_H

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_exact_diagonalization_advanced/advanced_ed_structures/phi_operators.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_exact_diagonalization_advanced/advanced_ed_structures/psi_state.h"

#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace DCA {
namespace ADVANCED_EXACT_DIAGONALIZATION {
// DCA::ADVANCED_EXACT_DIAGONALIZATION::

template <typename parameter_type, typename ed_options>  // N: size of bitset sequence
class symmetry_operation {
public:
  typedef typename ed_options::b_dmn b_dmn;
  typedef typename ed_options::s_dmn s_dmn;
  typedef typename ed_options::r_dmn r_dmn;
  typedef typename ed_options::k_dmn k_dmn;

public:
  symmetry_operation();

  void print();

  // Initialize with vector containing permutation of bits corresponding to symmetry operation
  void initialize(const std::vector<int>& permutation_vector);

  void execute(psi_state<parameter_type, ed_options>& Psi, bool change_coeffs = true);

  int get_order() {
    return order;
  }

  // For debugging
  std::string get_name() {
    return "rotation_90";
  }

private:
  int order;
  std::vector<int> permutation;
};

template <typename parameter_type, typename ed_options>
symmetry_operation<parameter_type, ed_options>::symmetry_operation() {}

template <typename parameter_type, typename ed_options>
void symmetry_operation<parameter_type, ed_options>::print() {
  std::stringstream ss;

  ss << "\n\n\t symmetry \n\n";
  ss << "\t\t order : " << order << "\n";
  for (int j = 0; j < permutation.size(); ++j)
    ss << "\t" << j << "\t" << permutation[j] << "\n";
  ss << "\n\n";

  std::cout << ss.str();
}

template <typename parameter_type, typename ed_options>
void symmetry_operation<parameter_type, ed_options>::initialize(
    const std::vector<int>& permutation_vector) {
  assert(permutation_vector.size() == b_dmn::dmn_size() * s_dmn::dmn_size() * r_dmn::dmn_size());

  permutation = permutation_vector;

  // Calculate order of permutation
  order = 0;

  std::vector<int> test(permutation.size());
  std::vector<int> test_permuted(permutation.size());
  std::vector<int> tmp(permutation.size());

  for (int i = 0; i < test.size(); ++i) {
    test[i] = i;
    test_permuted[i] = i;
  }

  do {
    for (int j = 0; j < permutation.size(); ++j)
      tmp[permutation[j]] = test_permuted[j];

    test_permuted.swap(tmp);
    ++order;
  } while (test_permuted != test);
}

template <typename parameter_type, typename ed_options>
void symmetry_operation<parameter_type, ed_options>::execute(psi_state<parameter_type, ed_options>& Psi,
                                                             bool change_coeffs) {
  // int sign  = 1;

  for (int i = 0; i < Psi.size(); ++i) {
    typename psi_state<parameter_type, ed_options>::phi_type phi_tmp(0);

    int sign = 1;

    for (int j = permutation.size() - 1; j >= 0; --j) {
      if (Psi.get_phi(i)[j]) {
        operators<parameter_type, ed_options>::create_at(permutation[j], phi_tmp, sign);
      }
    }

    Psi.get_phi(i) = phi_tmp;

    if (change_coeffs) {
      Psi.get_alpha(i) *= sign;
    }
  }

  Psi.sort();
}

}  // ADVANCED_EXACT_DIAGONALIZATION
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ADVANCED_ED_STRUCTURES_PSI_SYMMETRY_OPERATOR_H
