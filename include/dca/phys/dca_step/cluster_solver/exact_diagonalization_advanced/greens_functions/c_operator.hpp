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
// Helper class for SpGreensFunction and TpGreensFunction.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_GREENS_FUNCTIONS_C_OPERATOR_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_GREENS_FUNCTIONS_C_OPERATOR_HPP

namespace dca {
namespace phys {
namespace solver {
namespace ed {
// dca::phys::solver::ed::

struct c_operator {
  int index;
  int bsr_ind;
  bool creation;
};

}  // ed
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_GREENS_FUNCTIONS_C_OPERATOR_HPP
