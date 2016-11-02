// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements feynman_expansion_order_domain.hpp.

#include "dca/phys/dca_step/cluster_solver/ctaux/domains/feynman_expansion_order_domain.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

std::vector<int>& Feynman_expansion_order_domain::initialize_elements() {
  static std::vector<int> v(get_size());

  for (int i = 0; i < get_size(); i++)
    v[i] = i;

  return v;
}

}  // ctaux
}  // solver
}  // phys
}  // dca
