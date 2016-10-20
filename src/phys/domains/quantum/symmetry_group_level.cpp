// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
