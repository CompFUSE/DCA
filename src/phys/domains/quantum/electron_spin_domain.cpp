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
#include <cassert>

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

void electron_spin_domain::initialize_elements(std::vector<e_spin_states_type>& elements) {
  assert(elements.size() > 0 && elements.size() < 3);
  for (int i_spin = 0; i_spin < elements.size(); ++i_spin) {
    elements.push_back(i_spin ? e_DN : e_UP);
  }
}

std::vector<e_spin_states_type> electron_spin_domain::initialize_elements() {
  if (elements_.empty()) {
    elements_.push_back(e_DN);
    elements_.push_back(e_UP);
  }
}

}  // namespace domains
}  // namespace phys
}  // namespace dca
