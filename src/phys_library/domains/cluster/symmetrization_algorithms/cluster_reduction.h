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

#ifndef PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIZATION_ALGORITHMS_CLUSTER_REDUCTION_H
#define PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIZATION_ALGORITHMS_CLUSTER_REDUCTION_H

#include "phys_library/domains/cluster/symmetrization_algorithms/search_maximal_symmetry_group.h"
#include "phys_library/domains/cluster/symmetrization_algorithms/set_symmetry_matrices.h"

template <class base_cluster_type, class point_group>
class cluster_reduction {
public:
  cluster_reduction();

  void execute();
};

template <class base_cluster_type, class point_group>
cluster_reduction<base_cluster_type, point_group>::cluster_reduction() {}

template <class base_cluster_type, class point_group>
void cluster_reduction<base_cluster_type, point_group>::execute() {
  // cout << __FUNCTION__ << endl;

  search_maximal_symmetry_group<base_cluster_type, point_group, UNIT_CELL>::execute();

  search_maximal_symmetry_group<base_cluster_type, point_group, SUPER_CELL>::execute();

  set_symmetry_matrices<base_cluster_type>::execute();

  //   set_symmetry_matrices<base_cluster_type>::print_on_shell();
  //   throw std::logic_error(__FUNCTION__);
}

#endif  // PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIZATION_ALGORITHMS_CLUSTER_REDUCTION_H
