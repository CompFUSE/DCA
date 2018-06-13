// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class manages the symmetry of the cluster.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_SYMMETRY_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_SYMMETRY_HPP

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/cluster_definitions.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain_family.hpp"
#include "dca/phys/domains/cluster/dual_cluster.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/point_group_symmetry_domain.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <typename scalar_type, int D, CLUSTER_NAMES N, CLUSTER_REPRESENTATION R, CLUSTER_SHAPE S>
class cluster_symmetry<cluster_domain<scalar_type, D, N, R, S>> {
  const static int DIMENSION = D;

  const static CLUSTER_NAMES NAME = N;
  const static CLUSTER_REPRESENTATION REPRESENTATION = R;
  const static CLUSTER_SHAPE SHAPE = S;

  const static CLUSTER_REPRESENTATION DUAL_REPRESENTATION =
      dual_cluster<REPRESENTATION>::REPRESENTATION;

public:
  typedef cluster_domain_family<scalar_type, D, N, S> cluster_family_type;

  typedef cluster_domain<scalar_type, D, N, REPRESENTATION, S> this_type;
  typedef cluster_domain<scalar_type, D, N, DUAL_REPRESENTATION, S> dual_type;

  typedef func::dmn_0<domains::electron_band_domain> b_dmn_t;
  typedef func::dmn_0<this_type> c_dmn_t;

  typedef func::dmn_0<domains::point_group_symmetry_domain<domains::UNIT_CELL, cluster_family_type>>
      sym_unit_cell_dmn_t;
  typedef func::dmn_0<domains::point_group_symmetry_domain<domains::SUPER_CELL, cluster_family_type>>
      sym_super_cell_dmn_t;

  typedef func::dmn_variadic<func::dmn_variadic<c_dmn_t, b_dmn_t>, sym_super_cell_dmn_t> symmetry_matrix_dmn_t;

  static func::function<std::pair<int, int>, symmetry_matrix_dmn_t>& get_symmetry_matrix() {
    static func::function<std::pair<int, int>, symmetry_matrix_dmn_t> symmetry_matrix(
        "symmetry_matrix_super_cell");
    return symmetry_matrix;
  }
};

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_CLUSTER_SYMMETRY_HPP
