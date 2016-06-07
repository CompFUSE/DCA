// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// Description

#ifndef PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_2D_2D_HEXAGONAL_H
#define PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_2D_2D_HEXAGONAL_H

#include "dca/util/type_list.hpp"
#include "phys_library/domains/cluster/symmetries/symmetry_operations/2D/Cn.h"
#include "phys_library/domains/cluster/symmetries/symmetry_operations/2D/Sn.h"

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

#endif  // PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_2D_2D_HEXAGONAL_H
