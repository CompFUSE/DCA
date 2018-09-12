// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements the conversion between ErrorComputationType and string.

#include "dca/phys/error_computation_type.hpp"

#include <stdexcept>

namespace dca {
namespace phys {
// dca::phys::

ErrorComputationType stringToErrorComputationType(const std::string& str) {
  if (str == "NONE")
    return ErrorComputationType::NONE;
  else if (str == "STANDARD_DEVIATION")
    return ErrorComputationType::STANDARD_DEVIATION;
  else if (str == "JACK_KNIFE")
    return ErrorComputationType::JACK_KNIFE;
  else
    throw(std::logic_error("Invalid error computation type."));
}

std::string toString(const ErrorComputationType type) {
  switch (type) {
    case ErrorComputationType::NONE:
      return "NONE";
    case ErrorComputationType::STANDARD_DEVIATION:
      return "STANDARD_DEVIATION";
    case ErrorComputationType::JACK_KNIFE:
      return "JACK_KNIFE";
    default:
      throw(std::logic_error("Invalid error computation type."));
  }
}

}  // phys
}  // dca
