// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// D4.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_3D_D4_3D_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_3D_D4_3D_HPP

#include "dca/phys/domains/cluster/symmetries/point_group.hpp"
#include "dca/phys/domains/cluster/symmetries/symmetry_operations/3d/Cn_3d.hpp"
#include "dca/phys/domains/cluster/symmetries/symmetry_operations/3d/P_3d.hpp"
#include "dca/phys/domains/cluster/symmetries/symmetry_operations/3d/Sn_3d.hpp"
#include "dca/util/type_list.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

class D4_3D : public point_group<3> {
public:
  const static int size = 9;

  typedef P_3D inversion_type;

  typedef Cn_3D<0, 0, 1, 1, 4> Cn_3D_Z_1_4_type;
  typedef Cn_3D<0, 0, 1, 2, 4> Cn_3D_Z_2_4_type;
  typedef Cn_3D<0, 0, 1, 3, 4> Cn_3D_Z_3_4_type;

  typedef Sn_3D<2, 0, 8> Sn_3D_0_8_type;
  typedef Sn_3D<2, 1, 8> Sn_3D_1_8_type;
  typedef Sn_3D<2, 2, 8> Sn_3D_2_8_type;
  typedef Sn_3D<2, 3, 8> Sn_3D_3_8_type;
  typedef Sn_3D<2, 4, 8> Sn_3D_4_8_type;

  typedef dca::util::Typelist<inversion_type, Cn_3D_Z_1_4_type, Cn_3D_Z_2_4_type, Cn_3D_Z_3_4_type,
                              Sn_3D_0_8_type, Sn_3D_1_8_type, Sn_3D_2_8_type, Sn_3D_3_8_type, Sn_3D_4_8_type>
      point_group_type_list;
};

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_3D_D4_3D_HPP
