// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Empty class template for the data/matrices that a Monte Carlo walker has to store.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_TEMPLATE_MC_WALKER_DATA_HPP
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_TEMPLATE_MC_WALKER_DATA_HPP

#include "comp_library/linalg/linalg_device_types.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_template/qmci_type.hpp"

namespace DCA {
namespace QMCI {
// DCA::QMCI::

template <QmciType qmci_type, LIN_ALG::device_type device_t, typename parameters_type>
class MC_walker_data {
public:
  MC_walker_data(parameters_type& parameters);
};

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_TEMPLATE_MC_WALKER_DATA_HPP
