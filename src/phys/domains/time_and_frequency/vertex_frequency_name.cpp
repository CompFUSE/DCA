// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
