// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Product group action.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_PRODUCT_GROUP_ACTION_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_PRODUCT_GROUP_ACTION_HPP

#include "dca/phys/domains/cluster/symmetries/symmetry_operations/matrix_product.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <typename term_1, typename term_2>
class product_group_action {
public:
  typedef typename term_1::base_type base_type;
  typedef product_group_action<term_1, term_2> this_type;

  typedef typename term_1::base_type base_type_1;
  typedef typename term_2::base_type base_type_2;

  const static double* matrix() {
    const static double* matrix = matrix_product<term_1, term_2, base_type_1, base_type_2>::matrix();
    return matrix;
  }
};

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_SYMMETRIES_SYMMETRY_OPERATIONS_PRODUCT_GROUP_ACTION_HPP
