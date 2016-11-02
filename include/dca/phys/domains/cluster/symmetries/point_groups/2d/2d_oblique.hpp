// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
