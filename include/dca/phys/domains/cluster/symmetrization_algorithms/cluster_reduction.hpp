// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class performs the cluster reduction.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIZATION_ALGORITHMS_CLUSTER_REDUCTION_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIZATION_ALGORITHMS_CLUSTER_REDUCTION_HPP

#include "dca/phys/domains/cluster/symmetrization_algorithms/search_maximal_symmetry_group.hpp"
#include "dca/phys/domains/cluster/symmetrization_algorithms/set_symmetry_matrices.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <class base_cluster_type, class point_group>
class cluster_reduction {
public:
  cluster_reduction() {}

  void execute();
};

template <class base_cluster_type, class point_group>
void cluster_reduction<base_cluster_type, point_group>::execute() {
  // cout << __FUNCTION__ << endl;

  search_maximal_symmetry_group<base_cluster_type, point_group, domains::UNIT_CELL>::execute();

  search_maximal_symmetry_group<base_cluster_type, point_group, domains::SUPER_CELL>::execute();

  set_symmetry_matrices<base_cluster_type>::execute();

  //   set_symmetry_matrices<base_cluster_type>::print_on_shell();
  //   throw std::logic_error(__FUNCTION__);
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIZATION_ALGORITHMS_CLUSTER_REDUCTION_HPP
