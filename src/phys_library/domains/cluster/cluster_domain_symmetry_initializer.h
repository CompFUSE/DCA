// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef PHYS_LIBRARY_DOMAINS_CLUSTER_CLUSTER_DOMAIN_SYMMETRY_INITIALIZER_H
#define PHYS_LIBRARY_DOMAINS_CLUSTER_CLUSTER_DOMAIN_SYMMETRY_INITIALIZER_H

#include "phys_library/domains/cluster/cluster_domain_symmetry.h"
#include "phys_library/domains/cluster/symmetrization_algorithms/cluster_reduction.h"

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
class cluster_domain_symmetry_initializer<dmn_0<cluster_type>, point_group_type> {
  typedef typename cluster_symmetry<cluster_type>::cluster_family_type cluster_family_type;

public:
  static void execute() {
    cluster_reduction<cluster_family_type, point_group_type> cluster_reduction_obj;
    cluster_reduction_obj.execute();
  }
};

#endif  // PHYS_LIBRARY_DOMAINS_CLUSTER_CLUSTER_DOMAIN_SYMMETRY_INITIALIZER_H
