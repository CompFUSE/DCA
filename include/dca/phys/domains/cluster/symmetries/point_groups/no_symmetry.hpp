// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// No symmetry.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_NO_SYMMETRY_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_NO_SYMMETRY_HPP

#include "dca/phys/domains/cluster/symmetries/symmetry_operations/identity_group_operation.hpp"
#include "dca/util/type_list.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <int DIMENSION>
struct no_symmetry {
  typedef dca::util::Typelist<identity_group_operation<DIMENSION>> point_group_type_list;
};

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_POINT_GROUPS_NO_SYMMETRY_HPP
