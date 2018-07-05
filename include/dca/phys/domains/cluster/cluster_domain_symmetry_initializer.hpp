// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class initializes the cluster symmetry.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_DOMAIN_SYMMETRY_INITIALIZER_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_DOMAIN_SYMMETRY_INITIALIZER_HPP

#include "dca/function/domains.hpp"
#include "dca/phys/domains/cluster/cluster_symmetry.hpp"
#include "dca/phys/domains/cluster/symmetrization_algorithms/cluster_reduction.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <typename cluster_type, typename point_group_type>
class cluster_domain_symmetry_initializer {
  typedef typename cluster_symmetry<cluster_type>::cluster_family_type cluster_family_type;

public:
  static void execute() {
    cluster_reduction<cluster_family_type, point_group_type> cluster_reduction_obj;
    cluster_reduction_obj.execute();
  }
};

template <typename cluster_type, typename point_group_type>
class cluster_domain_symmetry_initializer<func::dmn_0<cluster_type>, point_group_type> {
  typedef typename cluster_symmetry<cluster_type>::cluster_family_type cluster_family_type;

public:
  static void execute() {
    cluster_reduction<cluster_family_type, point_group_type> cluster_reduction_obj;
    cluster_reduction_obj.execute();
  }
};

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_DOMAIN_SYMMETRY_INITIALIZER_HPP
