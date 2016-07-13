// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_2D_2D_SQUARE_H
#define PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_2D_2D_SQUARE_H

#include "dca/util/type_list.hpp"
#include "phys_library/domains/cluster/symmetries/symmetry_operations/2D/Cn.h"
#include "phys_library/domains/cluster/symmetries/symmetry_operations/2D/Sn.h"
#include "phys_library/domains/cluster/symmetries/symmetry_operations/identity_operation.h"

// Group actions
typedef Cn_2D<1, 4> Cn_2D_1_4_type;
typedef Cn_2D<2, 4> Cn_2D_2_4_type;
typedef Cn_2D<3, 4> Cn_2D_3_4_type;
typedef Cn_2D<4, 4> Cn_2D_4_4_type;

typedef Sn_2D<0, 8> Sn_2D_0_8_type;
typedef Sn_2D<1, 8> Sn_2D_1_8_type;
typedef Sn_2D<2, 8> Sn_2D_2_8_type;
typedef Sn_2D<3, 8> Sn_2D_3_8_type;
typedef Sn_2D<4, 8> Sn_2D_4_8_type;
typedef Sn_2D<5, 8> Sn_2D_5_8_type;
typedef Sn_2D<6, 8> Sn_2D_6_8_type;
typedef Sn_2D<7, 8> Sn_2D_7_8_type;
typedef Sn_2D<8, 8> Sn_2D_8_8_type;

// Point group: set of group actions
struct C4 {
  typedef dca::util::Typelist<Cn_2D_1_4_type, Cn_2D_2_4_type, Cn_2D_3_4_type, Cn_2D_4_4_type>
      point_group_type_list;
};

struct S4 {
  typedef dca::util::Typelist<Sn_2D_0_8_type, Sn_2D_2_8_type, Sn_2D_4_8_type, Sn_2D_6_8_type>
      point_group_type_list;
};

struct S4_plus {
  typedef dca::util::Typelist<identity_group_operation<2>, Cn_2D_2_4_type, Sn_2D_0_8_type, Sn_2D_2_8_type>
      point_group_type_list;
};

struct S8 {
  typedef dca::util::Typelist<Sn_2D_0_8_type, Sn_2D_1_8_type, Sn_2D_2_8_type, Sn_2D_3_8_type,
                              Sn_2D_4_8_type, Sn_2D_5_8_type, Sn_2D_6_8_type, Sn_2D_7_8_type,
                              identity_group_operation<2>>
      point_group_type_list;
};

struct D4 {
  typedef dca::util::Typelist<Cn_2D_1_4_type, Cn_2D_2_4_type, Cn_2D_3_4_type, Sn_2D_0_8_type,
                              Sn_2D_1_8_type, Sn_2D_2_8_type, Sn_2D_3_8_type, Sn_2D_4_8_type, Sn_2D_5_8_type,
                              Sn_2D_6_8_type, Sn_2D_7_8_type, identity_group_operation<2>>
      point_group_type_list;
};

#endif  // PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_2D_2D_SQUARE_H
