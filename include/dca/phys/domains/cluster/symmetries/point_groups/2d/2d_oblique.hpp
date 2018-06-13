// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// 2D oblique.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_2D_2D_OBLIQUE_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_2D_2D_OBLIQUE_HPP

#include "dca/phys/domains/cluster/symmetries/symmetry_operations/2d/Cn_2d.hpp"
#include "dca/util/type_list.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

// Group actions
typedef Cn_2D<1, 2> Cn_2D_1_2_type;

// Point group: set of group actions
typedef dca::util::Typelist<Cn_2D_1_2_type> C2;

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_2D_2D_OBLIQUE_HPP
