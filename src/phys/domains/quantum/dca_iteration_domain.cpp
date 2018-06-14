// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements dca_iteration_domain.hpp.

#include "dca/phys/domains/quantum/dca_iteration_domain.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

std::vector<int>& DCA_iteration_domain::initialize_elements() {
  static std::vector<int> v(get_size());

  for (int i = 0; i < get_size(); i++)
    v[i] = i;

  return v;
}

}  // domains
}  // phys
}  // dca
