// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Empty class template for a Monte Carlo accumulator.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_TEMPLATE_MC_ACCUMULATOR_HPP
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_TEMPLATE_MC_ACCUMULATOR_HPP

#include "comp_library/linalg/linalg_device_types.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_template/qmci_type.hpp"

namespace DCA {
namespace QMCI {
// DCA::QMCI::

template <QmciType qmci_type, LIN_ALG::device_type device_t, class parameters_type, class MOMS_type>
class MC_accumulator {
public:
  MC_accumulator(parameters_type& parameters_ref, MOMS_type& MOMS_ref);

  template <typename dca_info_struct_t>
  void finalize(dca_info_struct_t& dca_info_struct);

  void initialize();

  void accumulate();

  template <class stream_type>
  void to_JSON(stream_type& ss);
};

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_TEMPLATE_MC_ACCUMULATOR_HPP
