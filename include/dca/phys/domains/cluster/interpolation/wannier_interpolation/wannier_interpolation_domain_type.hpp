// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file defines the Wannier interpolation domain type.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_INTERPOLATION_WANNIER_INTERPOLATION_WANNIER_INTERPOLATION_DOMAIN_TYPE_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_INTERPOLATION_WANNIER_INTERPOLATION_WANNIER_INTERPOLATION_DOMAIN_TYPE_HPP

#include <dca/util/type_utils.hpp>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <typename dmn_type, typename type_input, typename type_output>
struct wannier_interpolation_domain_type {
  typedef typename dmn_type::this_type dmn_type_list;
  typedef typename dca::util::Swap<dmn_type_list, type_input, type_output>::type Result;
};

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_INTERPOLATION_WANNIER_INTERPOLATION_WANNIER_INTERPOLATION_DOMAIN_TYPE_HPP
