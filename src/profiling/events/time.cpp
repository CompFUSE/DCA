// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements time.hpp.

#include "dca/profiling/events/time.hpp"

namespace dca {
namespace profiling {
// dca::profiling::

Duration operator-(const Time& time, const Time& earlierTime) {
  return Duration(time, earlierTime);
}

}  // profiling
}  // dca
