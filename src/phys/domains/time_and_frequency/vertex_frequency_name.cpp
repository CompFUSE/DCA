// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements vertex_frequency_name.hpp.

#include "dca/phys/domains/time_and_frequency/vertex_frequency_name.hpp"
#include <stdexcept>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

std::string to_str(VERTEX_FREQUENCY_NAME NAME) {
  switch (NAME) {
    case COMPACT:
      return "COMPACT";
      break;

    case EXTENDED:
      return "EXTENDED";
      break;

    case COMPACT_POSITIVE:
      return "COMPACT_POSITIVE";
      break;

    case EXTENDED_POSITIVE:
      return "EXTENDED_POSITIVE";
      break;

    case EXTENDED_BOSONIC:
      return "EXTENDED_BOSONIC";
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }

  return "NONE";
}

}  // domains
}  // phys
}  // dca
