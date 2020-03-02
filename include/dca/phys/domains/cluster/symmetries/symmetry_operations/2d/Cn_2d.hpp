;// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Giovanni Balduzzi (gbalduzz@phys.itp.ethz.ch)
//
// Rotation over a 2*pi*n/m.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_2D_CN_2D_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_2D_CN_2D_HPP

#include "dca/phys/domains/cluster/symmetries/symmetry_operations/group_action.hpp"

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
    static std::array<double, 4> matrix;

    constexpr double c = std::cos(2 * M_PI * double(n) / double(m));
    constexpr double s = std::sin(2 * M_PI * double(n) / double(m));

    // in column major order !
    matrix[0] = c;
    matrix[1] = s;
    matrix[2] = -s;
    matrix[3] = c;

    return matrix.data();
  }
};

}  // namespace domains
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_2D_CN_2D_HPP
