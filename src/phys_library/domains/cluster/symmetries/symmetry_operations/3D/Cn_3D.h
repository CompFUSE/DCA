// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// Elements of the cyclic group C_n(m).

#ifndef PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_CN_3D_H
#define PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_CN_3D_H

#include "phys_library/domains/cluster/symmetries/symmetry_operations/group_action.h"
#include "phys_library/domains/cluster/symmetries/symmetry_operations/trigoniometric_ops/trig_ops.h"

template <int ux, int uy, int uz, int n, int m>
class Cn_3D : public group_action<3> {
public:
  typedef group_action<3> base_type;
  typedef Cn_3D<ux, uy, uz, n, m> this_type;

  Cn_3D(){};

  static double* matrix() {
    // rotation around an axis {ux,uy,uz} with angle th = 2*pi*n/m;

    double c = COSINE_EVAL<n, m>::value();  // cos(2*M_PI*double(n)/double(m));
    double s = SINE_EVAL<n, m>::value();    // sin(2*M_PI*double(n)/double(m));

    double Ux = double(ux) / sqrt(double(ux * ux + uy * uy + uz * uz));
    double Uy = double(uy) / sqrt(double(ux * ux + uy * uy + uz * uz));
    double Uz = double(uz) / sqrt(double(ux * ux + uy * uy + uz * uz));

    static double matrix[3 * 3] = {
        Ux * Ux + (1 - Ux * Ux) * c, Ux * Uy * (1 - c) + Uz * s,  Ux * Uz * (1 - c) - Uy * s,
        Ux * Uy * (1 - c) - Uz * s,  Uy * Uy + (1 - Uy * Uy) * c, Uy * Uz * (1 - c) + Ux * s,
        Ux * Uz * (1 - c) + Uy * s,  Uy * Uz * (1 - c) - Ux * s,  Uz * Uz + (1 - Uz * Uz) * c};
    return matrix;
  }

private:
};

#endif  // PHYS_LIBRARY_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_CN_3D_H
