// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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

enum ClusterSolverName { CT_AUX, SS_CT_HYB, ED_ADVANCED, HIGH_TEMPERATURE_SERIES };

}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_NAME_HPP
