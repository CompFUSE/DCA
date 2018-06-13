// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// 2D hexagonal.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_2D_2D_HEXAGONAL_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_2D_2D_HEXAGONAL_HPP

#include "dca/phys/domains/cluster/symmetries/symmetry_operations/2d/Cn_2d.hpp"
#include "dca/phys/domains/cluster/symmetries/symmetry_operations/2d/Sn_2d.hpp"
#include "dca/util/type_list.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

// Group actions
typedef Cn_2D<1, 3> Cn_2D_1_3_type;
typedef Cn_2D<2, 3> Cn_2D_2_3_type;

typedef Sn_2D<0, 3> Sn_2D_0_3_type;
typedef Sn_2D<1, 3> Sn_2D_1_3_type;
typedef Sn_2D<2, 3> Sn_2D_2_3_type;

typedef Cn_2D<1, 6> Cn_2D_1_6_type;
typedef Cn_2D<2, 6> Cn_2D_2_6_type;
typedef Cn_2D<3, 6> Cn_2D_3_6_type;
typedef Cn_2D<4, 6> Cn_2D_4_6_type;
typedef Cn_2D<5, 6> Cn_2D_5_6_type;

typedef Sn_2D<0, 6> Sn_2D_0_6_type;
typedef Sn_2D<1, 6> Sn_2D_1_6_type;
typedef Sn_2D<2, 6> Sn_2D_2_6_type;
typedef Sn_2D<3, 6> Sn_2D_3_6_type;
typedef Sn_2D<4, 6> Sn_2D_4_6_type;
typedef Sn_2D<5, 6> Sn_2D_5_6_type;

// Point group: set of group actions
struct C3 {
  typedef dca::util::Typelist<Cn_2D_1_3_type, Cn_2D_2_3_type> point_group_type_list;
};

struct S3 {
  typedef dca::util::Typelist<Sn_2D_0_3_type, Sn_2D_1_3_type, Sn_2D_2_3_type> point_group_type_list;
};

struct D3 {
  typedef dca::util::Typelist<Cn_2D_1_3_type, Cn_2D_2_3_type, Sn_2D_0_3_type, Sn_2D_1_3_type, Sn_2D_2_3_type>
      point_group_type_list;
};

struct C6 {
  typedef dca::util::Typelist<Cn_2D_1_6_type, Cn_2D_2_6_type, Cn_2D_3_6_type, Cn_2D_4_6_type, Cn_2D_5_6_type>
      point_group_type_list;
};

struct S6 {
  typedef dca::util::Typelist<Sn_2D_0_6_type, Sn_2D_1_6_type, Sn_2D_2_6_type, Sn_2D_3_6_type,
                              Sn_2D_4_6_type, Sn_2D_5_6_type>
      point_group_type_list;
};

struct D6 {
  typedef dca::util::Typelist<Cn_2D_1_6_type, Cn_2D_2_6_type, Cn_2D_3_6_type, Cn_2D_4_6_type,
                              Cn_2D_5_6_type, Sn_2D_0_6_type, Sn_2D_1_6_type, Sn_2D_2_6_type,
                              Sn_2D_3_6_type, Sn_2D_4_6_type, Sn_2D_5_6_type>
      point_group_type_list;
};

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_2D_2D_HEXAGONAL_HPP
