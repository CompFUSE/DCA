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

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ADVANCED_ED_STRUCTURES_PHI_OPERATIONS_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ADVANCED_ED_STRUCTURES_PHI_OPERATIONS_H

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_exact_diagonalization_advanced/advanced_ed_structures/phi_names.hpp"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_exact_diagonalization_advanced/advanced_ed_structures/phi_state.h"

namespace DCA {
namespace ADVANCED_EXACT_DIAGONALIZATION {
// DCA::ADVANCED_EXACT_DIAGONALIZATION::

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

}  // ADVANCED_EXACT_DIAGONALIZATION
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_ADVANCED_ED_STRUCTURES_PHI_OPERATIONS_H
