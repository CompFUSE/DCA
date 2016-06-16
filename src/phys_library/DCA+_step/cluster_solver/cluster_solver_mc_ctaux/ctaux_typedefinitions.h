// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// This class defines common types for the CT-AUX Monte Carlo integrator.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_TYPEDEFINITIONS_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_TYPEDEFINITIONS_H

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_template/mc_type_definitions.hpp"

namespace DCA {
namespace QMCI {
// DCA::QMCI::

template <class parameters_type, class MOMS_type>
class MC_type_definitions<CT_AUX_SOLVER, parameters_type, MOMS_type> {
public:
  // Types that define the profiling.
  typedef typename parameters_type::concurrency_type concurrency_type;
  typedef typename parameters_type::profiler_type profiler_type;
};

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_TYPEDEFINITIONS_H
