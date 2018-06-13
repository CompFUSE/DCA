// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Rotation over a 2*pi*n/m.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_2D_CN_2D_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_2D_CN_2D_HPP

#include "dca/phys/domains/cluster/symmetries/symmetry_operations/group_action.hpp"
#include "dca/phys/domains/cluster/symmetries/symmetry_operations/trigoniometric_ops/trig_ops.h"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

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

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_2D_CN_2D_HPP
