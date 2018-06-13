// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// 2D rectangular.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_2D_2D_RECTANGULAR_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_2D_2D_RECTANGULAR_HPP

#include "dca/phys/domains/cluster/symmetries/symmetry_operations/2d/Sn_2d.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_oblique.hpp"  // for C2
#include "dca/util/type_list.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

// Group actions
typedef Sn_2D<0, 4> Sn_2D_0_4_type;
typedef Sn_2D<1, 4> Sn_2D_1_4_type;

// Point group: set of group actions
typedef dca::util::Typelist<Sn_2D_0_4_type> Sx_2d;
typedef dca::util::Typelist<Sn_2D_1_4_type> Sy_2d;

typedef dca::util::Append<C2, dca::util::Append<Sx_2d, Sy_2d>::type>::type D2;

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_2D_2D_RECTANGULAR_HPP
