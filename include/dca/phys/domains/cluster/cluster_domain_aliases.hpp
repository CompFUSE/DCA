// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This class hold common type aliases related to cluster domains
//

#ifndef DCA_PHYS_CLUSTER_DOMAIN_ALIASES_HPP
#define DCA_PHYS_CLUSTER_DOMAIN_ALIASES_HPP

#include "dca/function/domains/dmn_0.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain_family.hpp"

namespace dca {
namespace phys {
// dca::phys::

template<int DIMENSION>
class ClusterDomainAliases {
public:
  // DCA cluster domains
  using RClusterType = domains::cluster_domain<double, DIMENSION, domains::CLUSTER,
					       domains::REAL_SPACE, domains::BRILLOUIN_ZONE>;
  using RClusterDmn =
      func::dmn_0<domains::cluster_domain<double, DIMENSION, domains::CLUSTER,
                                          domains::REAL_SPACE, domains::BRILLOUIN_ZONE>>;
  using KClusterType = domains::cluster_domain<double, DIMENSION, domains::CLUSTER,
					      domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>;
  using KClusterDmn =
      func::dmn_0<domains::cluster_domain<double, DIMENSION, domains::CLUSTER,
                                          domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;
  // Host cluster domains
  using RSpHostDmn =
      func::dmn_0<domains::cluster_domain<double, DIMENSION, domains::LATTICE_SP,
                                          domains::REAL_SPACE, domains::BRILLOUIN_ZONE>>;
  using KSpHostDmn =
      func::dmn_0<domains::cluster_domain<double, DIMENSION, domains::LATTICE_SP,
                                          domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;
  // Host vertex cluster domains
  using RTpHostDmn =
      func::dmn_0<domains::cluster_domain<double, DIMENSION, domains::LATTICE_TP,
                                          domains::REAL_SPACE, domains::BRILLOUIN_ZONE>>;
  using KTpHostDmn =
      func::dmn_0<domains::cluster_domain<double, DIMENSION, domains::LATTICE_TP,
                                          domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;
  using DcaClusterFamily =
      domains::cluster_domain_family<double, DIMENSION, domains::CLUSTER,
                                     domains::BRILLOUIN_ZONE>;
  using HostSpClusterFamily =
      domains::cluster_domain_family<double, DIMENSION, domains::LATTICE_SP,
                                     domains::BRILLOUIN_ZONE>;
  using HostTpClusterFamily =
      domains::cluster_domain_family<double, DIMENSION, domains::LATTICE_TP,
                                     domains::BRILLOUIN_ZONE>;
};

}  // phys
}  // dca

#endif  // DCA_PHYS_CLUSTER_DOMAIN_ALIASES_HPP
