// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides symmetry operations for psi states.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_BASIS_STATES_SYMMETRY_OPERATION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_BASIS_STATES_SYMMETRY_OPERATION_HPP

#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/basis_states/phi_fermionic_operators.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/basis_states/psi_state.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ed {
// dca::phys::solver::ed::

template <typename parameters_type, typename ed_options>  // N: size of bitset sequence
class symmetry_operation {
public:
  typedef typename ed_options::b_dmn b_dmn;
  typedef typename ed_options::s_dmn s_dmn;

  using CDA = ClusterDomainAliases<parameters_type::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using KClusterDmn = typename CDA::KClusterDmn;

public:
  symmetry_operation();

  void print();

  // Initialize with vector containing permutation of bits corresponding to symmetry operation
  void initialize(const std::vector<int>& permutation_vector);

  void execute(psi_state<parameters_type, ed_options>& Psi, bool change_coeffs = true);

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

template <typename parameters_type, typename ed_options>
symmetry_operation<parameters_type, ed_options>::symmetry_operation() {}

template <typename parameters_type, typename ed_options>
void symmetry_operation<parameters_type, ed_options>::print() {
  std::stringstream ss;

  ss << "\n\n\t symmetry \n\n";
  ss << "\t\t order : " << order << "\n";
  for (int j = 0; j < permutation.size(); ++j)
    ss << "\t" << j << "\t" << permutation[j] << "\n";
  ss << "\n\n";

  std::cout << ss.str();
}

template <typename parameters_type, typename ed_options>
void symmetry_operation<parameters_type, ed_options>::initialize(
    const std::vector<int>& permutation_vector) {
  assert(permutation_vector.size() == b_dmn::dmn_size() * s_dmn::dmn_size() * RClusterDmn::dmn_size());

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

template <typename parameters_type, typename ed_options>
void symmetry_operation<parameters_type, ed_options>::execute(psi_state<parameters_type, ed_options>& Psi,
                                                             bool change_coeffs) {
  // int sign  = 1;

  for (int i = 0; i < Psi.size(); ++i) {
    typename psi_state<parameters_type, ed_options>::phi_type phi_tmp(0);

    int sign = 1;

    for (int j = permutation.size() - 1; j >= 0; --j) {
      if (Psi.get_phi(i)[j]) {
        PhiFermionicOperators<parameters_type, ed_options>::create_at(permutation[j], phi_tmp, sign);
      }
    }

    Psi.get_phi(i) = phi_tmp;

    if (change_coeffs) {
      Psi.get_alpha(i) *= sign;
    }
  }

  Psi.sort();
}

}  // ed
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_BASIS_STATES_SYMMETRY_OPERATION_HPP
