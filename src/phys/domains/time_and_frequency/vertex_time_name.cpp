// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements vertex_time_name.hpp.

#include "dca/phys/domains/time_and_frequency/vertex_time_name.hpp"
#include <stdexcept>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

std::string to_str(VERTEX_TIME_NAME NAME) {
  switch (NAME) {
    case SP_TIME_DOMAIN:
      return "SP_TIME_DOMAIN";
      break;

    case SP_TIME_DOMAIN_POSITIVE:
      return "SP_TIME_DOMAIN_POSITIVE";
      break;

    case SP_TIME_DOMAIN_LEFT_ORIENTED:
      return "SP_TIME_DOMAIN_LEFT_ORIENTED";
      break;

    case SP_TIME_DOMAIN_LEFT_ORIENTED_POSITIVE:
      return "SP_TIME_DOMAIN_LEFT_ORIENTED_POSITIVE";
      break;

    case TP_TIME_DOMAIN:
      return "TP_TIME_DOMAIN";
      break;

    case TP_TIME_DOMAIN_POSITIVE:
      return "TP_TIME_DOMAIN_POSITIVE";
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }

  return "NONE";
}

}  // domains
}  // phys
}  // dca
