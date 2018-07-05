// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements symmetry_group_level.hpp.

#include "dca/phys/domains/quantum/symmetry_group_level.hpp"
#include <stdexcept>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

std::string to_str(symmetry_group_level name) {
  switch (name) {
    case UNIT_CELL:
      return "UNIT_CELL";
      break;

    case SUPER_CELL:
      return "SUPER_CELL";
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

}  // domains
}  // phys
}  // dca
