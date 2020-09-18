// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Method for extracting a single G0 sector.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_INTERPOLATION_SHRINK_G0_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_INTERPOLATION_SHRINK_G0_HPP

#include "dca/function/function.hpp"
#include "dca/function/domains/dmn_0.hpp"
#include "dca/function/domains/dmn_variadic.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace details {
// dca::phys::solver::details::

using TDmn = func::dmn_0<domains::time_domain>;
using SDmn = func::dmn_0<domains::electron_spin_domain>;
using BDmn = func::dmn_0<domains::electron_band_domain>;
using SDmn = func::dmn_0<domains::electron_spin_domain>;
using Nu = func::dmn_variadic<BDmn, SDmn>;

template <int dimension>
using RDmn = func::dmn_0<domains::cluster_domain<double, dimension, domains::CLUSTER,
                                                 domains::REAL_SPACE, domains::BRILLOUIN_ZONE>>;

template <int dimension>
using SpGreensFunction = func::function<double, func::dmn_variadic<Nu, Nu, RDmn<dimension>, TDmn>>;

template <int dimension>
auto shrinkG0(const SpGreensFunction<dimension>& G0) {
  func::function<double, func::dmn_variadic<BDmn, BDmn, RDmn<dimension>, TDmn>> g0_trimmed;
  const int s = 0;
  for (int b1 = 0; b1 < BDmn::dmn_size(); b1++)
    for (int b2 = 0; b2 < BDmn::dmn_size(); b2++)
      for (int r = 0; r < RDmn<dimension>::dmn_size(); r++)
        for (int t = 0; t < TDmn::dmn_size(); t++)
          g0_trimmed(b1, b2, r, t) = G0(b1, s, b2, s, r, t);
  return g0_trimmed;
}

}  // namespace details
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_INTERPOLATION_SHRINK_G0_HPP
