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

#ifndef PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_2D_2D_RECTANGULAR_H
#define PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_2D_2D_RECTANGULAR_H

#include "dca/util/type_list.hpp"
#include "phys_library/domains/cluster/symmetries/symmetry_operations/2D/Sn.h"
#include "phys_library/domains/cluster/symmetries/point_groups/2D/2D_oblique.h"  // for C2

// Group actions
typedef Sn_2D<0, 4> Sn_2D_0_4_type;
typedef Sn_2D<1, 4> Sn_2D_1_4_type;

// Point group: set of group actions
typedef dca::util::Typelist<Sn_2D_0_4_type> Sx_2d;
typedef dca::util::Typelist<Sn_2D_1_4_type> Sy_2d;

typedef dca::util::Append<C2, dca::util::Append<Sx_2d, Sy_2d>::type>::type D2;

#endif  // PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_2D_2D_RECTANGULAR_H
