// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Rotation around an axis {ux,uy,uz} with angle th = 2*pi*n/m.

#ifndef PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_P_3D_H
#define PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_P_3D_H

#include "phys_library/domains/cluster/symmetries/symmetry_operations/group_action.h"

class P_3D : public group_action<3> {
public:
  typedef group_action<3> base_type;
  typedef P_3D this_type;

  P_3D(){};

  const static double* matrix() {
    static double matrix[3 * 3] = {-1., 0., 0., 0., -1., 0., 0., 0., -1.};

    return matrix;
  }
};

#endif  // PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_P_3D_H
