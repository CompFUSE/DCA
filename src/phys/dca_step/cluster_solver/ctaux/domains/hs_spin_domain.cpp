// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements hs_spin_domain.hpp.

#include "dca/phys/dca_step/cluster_solver/ctaux/domains/hs_spin_domain.hpp"

#include <stdexcept>

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

std::vector<HS_spin_states_type> HS_spin_domain::initialize_elements() {
  static std::vector<HS_spin_states_type> v(0);

  v.push_back(HS_DN);
  v.push_back(HS_ZERO);
  v.push_back(HS_UP);

  return v;
}

int HS_spin_domain::to_coordinate(element_type spin) {
  switch (spin) {
    case HS_DN:
      return 0;
      break;

    case HS_ZERO:
      return 1;
      break;

    case HS_UP:
      return 2;
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

}  // ctaux
}  // solver
}  // phys
}  // dca
