// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Apply symmetries.

#ifndef PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIZATION_ALGORITHMS_APPLY_SYMMETRIES_H
#define PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIZATION_ALGORITHMS_APPLY_SYMMETRIES_H

#include "dca/util/type_list.hpp"
#include "phys_library/domains/cluster/symmetrization_algorithms/apply_symmetry.h"

template <class cluster_type, class point_group>
struct apply_symmetries {
  typedef typename point_group::point_group_type_list point_group_type_list;

  const static int INDEX = dca::util::Length<point_group_type_list>::value;

  template <class cluster_reduction_type>
  static void execute(cluster_reduction_type& cluster_reduction) {
    reduce_cluster_coordinates(cluster_reduction);
    reduce_cluster_pairs(cluster_reduction);
  }

  template <class cluster_reduction_type>
  static void reduce_cluster_coordinates(cluster_reduction_type& cluster_reduction) {
    apply_symmetry<cluster_type, point_group_type_list, INDEX - 1>::reduce_cluster_coordinates(
        cluster_reduction);
  }

  template <class cluster_reduction_type>
  static void reduce_cluster_pairs(cluster_reduction_type& cluster_reduction) {
    for (int i = cluster_type::get_size() - 1; i >= 0; i--) {
      if (i == cluster_reduction.get_irreducible_index(i))  // is an irreducible point
      {
        for (int j = 0; j < cluster_type::get_size(); j++) {
          cluster_reduction.get_irreducible_pair(i, j).first = i;
          cluster_reduction.get_irreducible_pair(i, j).second = j;
        }
      }
      else {
        apply_symmetry<cluster_type, point_group_type_list, INDEX - 1>::reduce_cluster_pair(
            i, cluster_reduction);
      }
    }
  }
};

template <class cluster_type>
struct apply_symmetries<cluster_type, PsiMag_symmetry_2D> {
  template <class cluster_reduction_type>
  static void execute(cluster_reduction_type& cluster_reduction) {}
};
#endif  // PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIZATION_ALGORITHMS_APPLY_SYMMETRIES_H
