// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class searches for the maximal symmetry group.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIZATION_ALGORITHMS_SEARCH_MAXIMAL_SYMMETRY_GROUP_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIZATION_ALGORITHMS_SEARCH_MAXIMAL_SYMMETRY_GROUP_HPP

#include "dca/phys/domains/cluster/symmetrization_algorithms/search_symmetry_group.hpp"
#include "dca/phys/domains/quantum/point_group_symmetry_domain.hpp"
#include "dca/phys/domains/quantum/symmetry_group_level.hpp"
#include "dca/util/type_list.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <class base_cluster_type, class point_group, domains::symmetry_group_level_type symmetry_group_level>
class search_maximal_symmetry_group {
public:
  typedef typename point_group::point_group_type_list point_group_type_list;

  const static int DIMENSION = base_cluster_type::DIMENSION;
  const static int MAX_SIZE = dca::util::Length<point_group_type_list>::value;

  //   typedef r_cluster<FULL, base_cluster_type> r_cluster_type;
  //   typedef k_cluster<FULL, base_cluster_type> k_cluster_type;

  typedef domains::point_group_symmetry_domain<symmetry_group_level, base_cluster_type> sym_dmn_t;

public:
  static void execute() {
    sym_dmn_t::DIMENSION = DIMENSION;
    sym_dmn_t::get_size() = 0;

    search_symmetry_group<base_cluster_type, point_group_type_list, symmetry_group_level,
                          MAX_SIZE - 1>::execute();

    //     sym_dmn_t::to_JSON(std::cout);
    //     throw std::logic_error(__FUNCTION__);
  }
};

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIZATION_ALGORITHMS_SEARCH_MAXIMAL_SYMMETRY_GROUP_HPP
