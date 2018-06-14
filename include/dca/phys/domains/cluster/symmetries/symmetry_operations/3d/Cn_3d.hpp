// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Elements of the cyclic group C_n(m).

#ifndef DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_CN_3D_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_CN_3D_HPP

#include "dca/phys/domains/cluster/symmetries/symmetry_operations/group_action.hpp"
#include "dca/phys/domains/cluster/symmetries/symmetry_operations/trigoniometric_ops/trig_ops.h"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

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
};

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_CN_3D_HPP
