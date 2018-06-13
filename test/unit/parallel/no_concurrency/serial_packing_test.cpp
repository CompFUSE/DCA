// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests serial_packing.hpp.

#include "dca/parallel/no_concurrency/serial_packing.hpp"

#include <cstdint>

#include "gtest/gtest.h"

#include "dca/function/domains/dmn.hpp"
#include "dca/function/domains/dmn_0.hpp"

TEST(SerialPackingTest, All) {
  dca::parallel::SerialPacking packing_interface;

  // Fundamental type
  int32_t i = 314;
  EXPECT_EQ(4, packing_interface.get_buffer_size(i));

  // std::string
  std::string s = "hello";
  EXPECT_EQ(5, packing_interface.get_buffer_size(s));

  // std::vector
  std::vector<int32_t> v{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  EXPECT_EQ(40, packing_interface.get_buffer_size(v));  // 10*4 = 40

  // function
  using Domain = dca::func::dmn_0<dca::func::dmn<3, int>>;
  dca::func::function<int32_t, Domain> f;
  EXPECT_EQ(12, packing_interface.get_buffer_size(f));  // 3*4 = 12
}
