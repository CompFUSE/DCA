// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Reflection over an axis which has an angle 2*pi*n/m with the x-axis.

#ifndef PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_IDENTITY_OPERATION_H
#define PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_IDENTITY_OPERATION_H

#include "phys_library/domains/cluster/symmetries/symmetry_operations/group_action.h"

template <int DIMENSION>
class identity_group_operation : public group_action<DIMENSION> {};

template <>
class identity_group_operation<2> : public group_action<2> {
public:
  typedef group_action<2> base_type;
  typedef identity_group_operation<2> this_type;

  identity_group_operation(){};

  const static double* matrix() {
    const static double matrix[2 * 2] = {1., 0., 0., 1.};
    return matrix;
  }
};

template <>
class identity_group_operation<3> : public group_action<3> {
public:
  typedef group_action<3> base_type;
  typedef identity_group_operation<3> this_type;

  identity_group_operation(){};

  const static double* matrix() {
    const static double matrix[3 * 3] = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
    return matrix;
  }
};

#endif  // PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_IDENTITY_OPERATION_H
