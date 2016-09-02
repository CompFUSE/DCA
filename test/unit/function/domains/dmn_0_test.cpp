// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi(gbalduzz@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests dmn_0.h

#include "comp_library/function_library/domains/special_domains/dmn_0.h"
#include "gtest/gtest.h"
#include "comp_library/function_library/domains/special_domains/dmn.h"

TEST(Dmn0Test, GetName) {
  using Domain1 = dmn_0<dmn<3, unsigned int>>;
  EXPECT_EQ("dmn_0<dmn<3, unsigned int>>", Domain1::get_name());

  using Domain2 = dmn_0<dmn<42, double>>;
  EXPECT_EQ("dmn_0<dmn<42, double>>", Domain2::get_name());
}
