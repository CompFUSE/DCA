// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi(gbalduzz@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests dmn_0.hpp

#include "dca/function/domains/dmn_0.hpp"
#include "gtest/gtest.h"
#include "dca/function/domains/dmn.hpp"

TEST(Dmn0Test, GetName) {
  using Domain1 = dca::func::dmn_0<dca::func::dmn<3, unsigned int>>;
  EXPECT_EQ("dmn_0<dca::func::dmn<3, unsigned int>>", Domain1::get_name());

  using Domain2 = dca::func::dmn_0<dca::func::dmn<42, double>>;
  EXPECT_EQ("dmn_0<dca::func::dmn<42, double>>", Domain2::get_name());
}
