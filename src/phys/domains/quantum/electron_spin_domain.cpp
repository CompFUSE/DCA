// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
