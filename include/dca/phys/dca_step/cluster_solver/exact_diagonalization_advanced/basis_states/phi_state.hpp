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
// This file provides a single basis state of the occupation number basis.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_BASIS_STATES_PHI_STATE_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_BASIS_STATES_PHI_STATE_HPP

#include <vector>
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/basis_states/phi_names.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ed {
// dca::phys::solver::ed::

// Class template declaration
template <typename parameter_type, typename ed_options,
          phi_names phi_name>  // N: size of bitset sequence
struct phi_state {};

// Phi singlet specialization
template <typename parameter_type, typename ed_options>  // N: size of bitset sequence
struct phi_state<parameter_type, ed_options, PHI_SINGLET> {
public:
  typedef typename ed_options::scalar_type scalar_type;
  typedef typename ed_options::complex_type complex_type;

  typedef typename ed_options::phi_type phi_type;

public:
  phi_type phi;
  complex_type alpha;
};

// Phi multiplet specialization
template <typename parameter_type, typename ed_options>  // N: size of bitset sequence
struct phi_state<parameter_type, ed_options, PHI_MULTIPLET> {
public:
  typedef typename ed_options::scalar_type scalar_type;
  typedef typename ed_options::complex_type complex_type;

  typedef typename ed_options::phi_type phi_type;

public:
  void sort();

public:
  phi_type phi;

  std::vector<int> index;
  std::vector<complex_type> alpha;
};

template <typename parameter_type, typename ed_options>
void phi_state<parameter_type, ed_options, PHI_MULTIPLET>::sort() {
  std::vector<int> sorted_index;
  std::vector<complex_type> sorted_alpha;

  for (int i = 0; i < index.size(); ++i) {
    int idx = 0;

    while (idx < sorted_index.size() && index[i] >= sorted_index[idx]) {
      ++idx;
    }

    sorted_index.insert(sorted_index.begin() + idx, index[i]);
    sorted_alpha.insert(sorted_alpha.begin() + idx, alpha[i]);
  }

  index.swap(sorted_index);
  alpha.swap(sorted_alpha);
}

}  // ed
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_BASIS_STATES_PHI_STATE_HPP
