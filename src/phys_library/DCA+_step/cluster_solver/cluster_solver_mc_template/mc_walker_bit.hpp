// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Empty class template for the built-in-test (BIT) of a Monte Carlo walker.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_TEMPLATE_MC_WALKER_BIT_HPP
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_TEMPLATE_MC_WALKER_BIT_HPP

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_template/qmci_type.hpp"

namespace DCA {
namespace QMCI {
// DCA::QMCI::

template <QmciType qmci_type, class parameters_type, class MOMS_type>
class MC_walker_BIT {
  // typedef MultiOrbitalMultiSiteStructure<parameters_type> MOMS_type;

public:
  MC_walker_BIT(parameters_type& parameters_ref, MOMS_type& MOMS_ref);
};

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_TEMPLATE_MC_WALKER_BIT_HPP
