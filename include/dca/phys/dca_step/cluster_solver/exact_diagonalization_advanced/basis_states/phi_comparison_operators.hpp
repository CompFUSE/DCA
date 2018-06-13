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
// This file provides comparison operators for phi states.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_BASIS_STATES_PHI_COMPARISON_OPERATORS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_BASIS_STATES_PHI_COMPARISON_OPERATORS_HPP

#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/basis_states/phi_names.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/basis_states/phi_state.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ed {
// dca::phys::solver::ed::

template <typename parameter_type, typename ed_options, phi_names phi_name>
bool operator<(const phi_state<parameter_type, ed_options, phi_name>& phi_obj1,
               const phi_state<parameter_type, ed_options, phi_name>& phi_obj2) {
  return phi_obj1.phi.to_ulong() < phi_obj2.phi.to_ulong();
}

template <typename parameter_type, typename ed_options>
bool operator==(const phi_state<parameter_type, ed_options, PHI_SINGLET>& phi_obj1,
                const phi_state<parameter_type, ed_options, PHI_SINGLET>& phi_obj2) {
  return (phi_obj1.phi == phi_obj2.phi and
          abs(phi_obj1.alpha - phi_obj2.alpha) < ed_options::get_epsilon());
}

}  // ed
}  // solver
}  // phys
}  // solver

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_BASIS_STATES_PHI_COMPARISON_OPERATORS_HPP
