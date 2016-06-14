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

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ADVANCED_ED_STRUCTURES_PSI_OPERATIONS_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ADVANCED_ED_STRUCTURES_PSI_OPERATIONS_H

#include <algorithm>
#include <complex>

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_exact_diagonalization_advanced/advanced_ed_structures/psi_state.h"

namespace DCA {
namespace ADVANCED_EXACT_DIAGONALIZATION {
// DCA::ADVANCED_EXACT_DIAGONALIZATION::

template <typename scalar_type, typename parameter_type, typename ed_options>
psi_state<parameter_type, ed_options> operator*(const scalar_type d,
                                                const psi_state<parameter_type, ed_options>& psi) {
  psi_state<parameter_type, ed_options> psi_new = psi;

  for (int i = 0; i < psi.size(); ++i) {
    psi_new.get_alpha(i) *= d;
  }

  return psi_new;
}

template <typename parameter_type, typename ed_options>
bool operator<(const psi_state<parameter_type, ed_options>& psi_1,
               const psi_state<parameter_type, ed_options>& psi_2) {
  int size_1 = psi_1.size();
  int size_2 = psi_2.size();

  for (int i = 0; i < std::min(size_1, size_2); ++i) {
    if (psi_1.get_phi(i).to_ulong() < psi_2.get_phi(i).to_ulong())
      return true;
    else if (psi_1.get_phi(i).to_ulong() > psi_2.get_phi(i).to_ulong())
      return false;
  }

  if (size_1 < size_2)
    return true;
  else
    return false;
}

template <typename parameter_type, typename ed_options>
bool less_magnetization(const psi_state<parameter_type, ed_options>& psi_1,
                        const psi_state<parameter_type, ed_options>& psi_2) {
  return psi_1.magnetization() < psi_2.magnetization();
}

template <typename parameter_type, typename ed_options>
bool possible_hopping(const typename psi_state<parameter_type, ed_options>::phi_type& phi1,
                      const typename psi_state<parameter_type, ed_options>::phi_type& phi2) {
  int diff = (phi1 ^ phi2).count();
  return (diff == 2);
}

// up to a phase?
// psi_1 and psi_2 first must be sorted
template <typename parameter_type, typename ed_options>
bool operator==(const psi_state<parameter_type, ed_options>& psi_1,
                const psi_state<parameter_type, ed_options>& psi_2) {
  return (psi_1.get_phi_obj() == psi_2.get_phi_obj());
}

template <typename parameter_type, typename ed_options>
typename ed_options::complex_type scalar_product(const psi_state<parameter_type, ed_options>& psi_1,
                                                 const psi_state<parameter_type, ed_options>& psi_2) {
  typedef typename ed_options::scalar_type scalar_type;
  typedef typename ed_options::complex_type complex_type;

  complex_type res = 0.;

  for (int i1 = 0; i1 < psi_1.size(); ++i1) {
    for (int i2 = 0; i2 < psi_2.size(); ++i2) {
      res += std::conj(psi_1.get_alpha(i1)) * psi_2.get_alpha(i2) *
             scalar_type(psi_1.get_phi(i1) == psi_2.get_phi(i2));
    }
  }

  return res;
}

}  // ADVANCED_EXACT_DIAGONALIZATION
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ADVANCED_ED_STRUCTURES_PSI_OPERATIONS_H
