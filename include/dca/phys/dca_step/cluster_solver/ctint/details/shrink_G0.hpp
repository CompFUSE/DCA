// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Method for extracting a single G0 sector.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_DETAILS_SHRINK_G0_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_DETAILS_SHRINK_G0_HPP

#include "dca/function/function.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/domains/common_domains.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
namespace details {
// dca::phys::solver::ctint::details::

template <class Rdmn>
auto shrinkG0(const func::function<double, func::dmn_variadic<Nu, Nu, Rdmn, Tdmn>>& G0) {
  func::function<double, func::dmn_variadic<Bdmn, Bdmn, Rdmn, Tdmn>> g0_trimmed;
  const int s = 0;
  for (int b1 = 0; b1 < Bdmn::dmn_size(); b1++)
    for (int b2 = 0; b2 < Bdmn::dmn_size(); b2++)
      for (int r = 0; r < Rdmn::dmn_size(); r++)
        for (int t = 0; t < Tdmn::dmn_size(); t++)
          g0_trimmed(b1, b2, r, t) = G0(b1, s, b2, s, r, t);
  return g0_trimmed;
}

}  // namespace details
}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_DETAILS_SHRINK_G0_HPP
