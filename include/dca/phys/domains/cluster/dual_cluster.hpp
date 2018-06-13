// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class provides the dual cluster of the cluster passed as template argument.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_DUAL_CLUSTER_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_DUAL_CLUSTER_HPP

#include "dca/phys/domains/cluster/cluster_definitions.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <CLUSTER_REPRESENTATION R>
struct dual_cluster {};

template <>
struct dual_cluster<MOMENTUM_SPACE> {
  const static CLUSTER_REPRESENTATION REPRESENTATION = REAL_SPACE;
};

template <>
struct dual_cluster<REAL_SPACE> {
  const static CLUSTER_REPRESENTATION REPRESENTATION = MOMENTUM_SPACE;
};

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_DUAL_CLUSTER_HPP
