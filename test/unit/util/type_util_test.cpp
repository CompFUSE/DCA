// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests the type_util helper methods.

#include "dca/util/type_utils.hpp"

#include "gtest/gtest.h"

TEST(TypeUtils, if_all){
  static_assert(not dca::util::if_all<1,1,0,1>::value, "Does not eval to false");
  static_assert(dca::util::if_all<1,1,1,1>::value, "Does not eval to true");
}
