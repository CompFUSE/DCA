// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements electron_spin_domain.hpp.

#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include <stdexcept>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

int electron_spin_domain::to_coordinate(element_type spin) {
  switch (spin) {
    case e_DN:
      return 0;
      break;

    case e_UP:
      return 1;
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

std::vector<e_spin_states_type> electron_spin_domain::initialize_elements() {
  static std::vector<e_spin_states_type> v(0);

  v.push_back(e_DN);
  v.push_back(e_UP);

  return v;
}

}  // domains
}  // phys
}  // dca
