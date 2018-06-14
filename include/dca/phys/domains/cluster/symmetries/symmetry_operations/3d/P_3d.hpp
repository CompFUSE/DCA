// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Rotation around an axis {ux,uy,uz} with angle th = 2*pi*n/m.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_P_3D_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_P_3D_HPP

#include "dca/phys/domains/cluster/symmetries/symmetry_operations/group_action.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

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

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_3D_P_3D_HPP
