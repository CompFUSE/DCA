// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Tyler Sax
//
// The two 2D crystallographic holohedries and the candidate pool they generate.
//
// D4 (square holohedry, order 8) and D6 (hexagonal holohedry, order 12) are the two maximal point
// groups in 2D. By the crystallographic restriction, every 2D crystallographic point group is a
// subgroup of one of them, so they are the only primitives needed: the lower-symmetry groups
// (rectangular D2, oblique C2, trigonal C3/D3, ...) are not declared per model but derived at
// runtime from the pool by the geometry filter (search_symmetry_group: is_lattice_compatible +
// is_unit_cell_compatible) and the H0-invariance gate.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_2D_HOLOHEDRIES_2D_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_2D_HOLOHEDRIES_2D_HPP

#include "dca/phys/domains/cluster/symmetries/symmetry_operations/2d/Cn_2d.hpp"
#include "dca/phys/domains/cluster/symmetries/symmetry_operations/2d/Sn_2d.hpp"
#include "dca/phys/domains/cluster/symmetries/symmetry_operations/identity_group_operation.hpp"
#include "dca/util/type_list.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

// Square holohedry (order 8): the three C4 rotations, four mirrors (reflection lines at n*pi/4,
// n={0...3}), and the identity.
struct D4 {
  using point_group_type_list =
      dca::util::Typelist<Sn_2D<0, 8>, Sn_2D<1, 8>, Sn_2D<2, 8>, Sn_2D<3, 8>, Cn_2D<1, 4>,
                          Cn_2D<2, 4>, Cn_2D<3, 4>, identity_group_operation<2>>;
};

// Hexagonal holohedry (order 12): the five C6 rotations, six mirrors (reflection lines at n*pi/6,
// n = 0..5), and the identity.
struct D6 {
  using point_group_type_list =
      dca::util::Typelist<Cn_2D<1, 6>, Cn_2D<2, 6>, Cn_2D<3, 6>, Cn_2D<4, 6>, Cn_2D<5, 6>,
                          Sn_2D<0, 12>, Sn_2D<1, 12>, Sn_2D<2, 12>, Sn_2D<3, 12>, Sn_2D<4, 12>,
                          Sn_2D<5, 12>, identity_group_operation<2>>;
};

// Universal 2D candidate pool = square holohedry D4  U  hexagonal holohedry D6. Overlapping ops
// (identity, the C2 rotation) appear twice; search_symmetry_group::is_duplicate collapses them at
// populate time, so the union need not be hand-deduplicated.
struct holohedry_pool_2D {
  using point_group_type_list =
      dca::util::Append<D4::point_group_type_list, D6::point_group_type_list>::type;
};

}  // namespace domains
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_2D_HOLOHEDRIES_2D_HPP
