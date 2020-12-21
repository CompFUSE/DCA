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
// Rotation over a 2*pi*n/m.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_2D_CN_2D_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_2D_CN_2D_HPP

#include <array>

#include "dca/phys/domains/cluster/symmetries/symmetry_operations/group_action.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <int n, int m>
class Cn_2D : public group_action<2> {
public:
  using base_type = group_action<2>;
  using this_type = Cn_2D<n, m>;

  Cn_2D(){};

  static const double* matrix() {
    static const auto matrix = init();
    return matrix.data();
  }

private:
  static std::array<double, 4> init() {
    std::array<double, 4> matrix;

    double c = std::cos(2 * M_PI * double(n) / double(m));
    double s = std::sin(2 * M_PI * double(n) / double(m));

    // in column major order !
    matrix[0] = c;
    matrix[1] = s;
    matrix[2] = -s;
    matrix[3] = c;

    return matrix;
  }
};

}  // namespace domains
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_2D_CN_2D_HPP
