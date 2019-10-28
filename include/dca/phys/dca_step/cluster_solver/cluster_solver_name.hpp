// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file defines the different types of cluster solvers.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_NAME_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_NAME_HPP

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

enum ClusterSolverName { CT_AUX, SS_CT_HYB, CT_INT, ED_ADVANCED, HIGH_TEMPERATURE_SERIES };

}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_NAME_HPP
