// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file defines the various types of cluster solvers.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_NAME_HPP
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_NAME_HPP

namespace DCA {

enum ClusterSolverName {
  HIGH_TEMPERATURE_SERIES,
  ED_CLUSTER_SOLVER,
  ADVANCED_ED_CLUSTER_SOLVER,
  CT_AUX_CLUSTER_SOLVER,
  SS_CT_HYB,
  ANALYSIS_INTERPOLATION
};

}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_NAME_HPP
