// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
//

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_DOMAINS_COMMON_DOMAINS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_DOMAINS_COMMON_DOMAINS_HPP

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
namespace ctint {
// dca::phys::solver::ctint::

using dca::func::dmn_0;
using dca::func::dmn_variadic;
using Tdmn = dmn_0<domains::time_domain>;
using Wdmn = dmn_0<domains::frequency_domain>;
using Sdmn = dmn_0<domains::electron_spin_domain>;
using Bdmn = dmn_0<domains::electron_band_domain>;
using Sdmn = dmn_0<domains::electron_spin_domain>;
using Nu = dmn_variadic<Bdmn, Sdmn>;


template<int DIMENSION>
using Rdmn_ = dmn_0<domains::cluster_domain<double, DIMENSION, domains::CLUSTER, domains::REAL_SPACE,
                                           domains::BRILLOUIN_ZONE>>;
template<int DIMENSION>
using Kdmn_ = dmn_0<domains::cluster_domain<double, DIMENSION, domains::CLUSTER,
                                           domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;


}  // ctint
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_DOMAINS_COMMON_DOMAINS_HPP
