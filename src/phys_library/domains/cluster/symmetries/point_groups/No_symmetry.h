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

#ifndef PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_NO_SYMMETRY_H
#define PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_NO_SYMMETRY_H

#include "dca/util/type_list.hpp"
#include "phys_library/domains/cluster/symmetries/symmetry_operations/identity_operation.h"

template <int DIMENSION>
struct no_symmetry {
  typedef dca::util::Typelist<identity_group_operation<DIMENSION>> point_group_type_list;
};

#endif  // PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_NO_SYMMETRY_H
