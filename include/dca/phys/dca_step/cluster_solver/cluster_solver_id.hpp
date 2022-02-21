// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// This file defines the different types of cluster solvers.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_NAME_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_NAME_HPP

namespace dca {
enum class ClusterSolverId { CT_AUX, SS_CT_HYB, CT_INT, ED_ADVANCED, HIGH_TEMPERATURE_SERIES };
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_NAME_HPP
