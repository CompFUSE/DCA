// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
