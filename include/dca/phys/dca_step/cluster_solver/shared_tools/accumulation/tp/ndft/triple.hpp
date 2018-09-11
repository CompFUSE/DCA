// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Small class used to sort the configuration.

#ifndef DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_TRIPLE_HPP
#define DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_TRIPLE_HPP

namespace dca {
namespace phys {
namespace solver {
namespace accumulator {
namespace details {
// dca::phys::solver::accumulator::details::

template <typename Real>
struct Triple {
  bool operator<(const Triple& rhs) const {
    return orbital == rhs.orbital ? idx < rhs.idx : orbital < rhs.orbital;
  }

  int orbital;
  int idx;
  Real tau;
};

}  // details
}  // accumulator
}  // solver
}  // phys
}  // dca

#endif  // DCA_INCLUDE_DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_ACCUMULATION_TP_TRIPLE_HPP
