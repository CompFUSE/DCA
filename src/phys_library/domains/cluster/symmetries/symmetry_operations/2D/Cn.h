// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Rotation over a 2*pi*n/m.

#ifndef PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_2D_CN_H
#define PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_2D_CN_H

#include "phys_library/domains/cluster/symmetries/symmetry_operations/group_action.h"
#include "phys_library/domains/cluster/symmetries/symmetry_operations/trigoniometric_ops/trig_ops.h"

template <int n, int m>
class Cn_2D : public group_action<2> {
public:
  typedef group_action<2> base_type;
  typedef Cn_2D<n, m> this_type;

  Cn_2D(){};

  static double* matrix() {
    static double* matrix = init();
    return matrix;
  }

private:
  static double* init() {
    static double* matrix = new double[4];

    double c = COSINE_EVAL<n, m>::value();  // cos(2.*M_PI*theta);
    double s = SINE_EVAL<n, m>::value();    // sin(2.*M_PI*theta);

    // in column major order !
    matrix[0] = c;
    matrix[1] = s;
    matrix[2] = -s;
    matrix[3] = c;

    return matrix;
  }
};

#endif  // PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_2D_CN_H
