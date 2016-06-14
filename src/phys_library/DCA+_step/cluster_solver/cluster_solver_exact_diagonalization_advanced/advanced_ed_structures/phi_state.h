// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// Description

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ADVANCED_ED_STRUCTURES_PHI_STATE_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ADVANCED_ED_STRUCTURES_PHI_STATE_H

#include <vector>
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_exact_diagonalization_advanced/advanced_ed_structures/phi_names.hpp"

namespace DCA {
namespace ADVANCED_EXACT_DIAGONALIZATION {
// DCA::ADVANCED_EXACT_DIAGONALIZATION::

template <typename parameter_type, typename ed_options,
          phi_names phi_name>  // N: size of bitset sequence
struct phi_state {};

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

}  // ADVANCED_EXACT_DIAGONALIZATION
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ADVANCED_ED_STRUCTURES_PHI_STATE_H
