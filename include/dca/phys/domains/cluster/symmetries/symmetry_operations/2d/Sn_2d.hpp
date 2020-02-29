// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Giovanni Balduzzi (gbalduzz@phys.itp.ethz.ch)
//
// Reflection over an axis which has an angle 2*pi*n/m with the x-axis.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_2D_SN_2D_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_2D_SN_2D_HPP

#include "dca/phys/domains/cluster/symmetries/symmetry_operations/group_action.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <int n, int m>
class Sn_2D : public group_action<2> {
public:
  typedef group_action<2> base_type;
  typedef Sn_2D<n, m> this_type;

  Sn_2D(){};

  static double* matrix() {
    static double* matrix = init();
    return matrix;
  }

private:
  static double* init() {
    static std::array<double, 4> matrix;

    constexpr double c = std::cos(2 * M_PI * double(n) / double(m));
    constexpr double s = std::sin(2 * M_PI * double(n) / double(m));

    // in column major order !
    matrix[0] = c * c - s * s;
    matrix[1] = 2 * s * c;
    matrix[2] = 2 * s * c;
    matrix[3] = -(c * c - s * s);

    return matrix.data();
  }
};

}  // namespace domains
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_2D_SN_2D_HPP
