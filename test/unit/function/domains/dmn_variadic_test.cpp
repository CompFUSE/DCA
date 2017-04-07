// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests dmn_variadic.hpp.

#include "dca/function/domains/dmn_variadic.hpp"
#include "gtest/gtest.h"
#include "dca/function/domains/dmn.hpp"
#include "dca/function/domains/dmn_0.hpp"

#ifndef NDEBUG

TEST(DmnVariadicTest, CheckIndices) {
  using dca::func::dmn_0;
  using dca::func::dmn;

  using SubDmn1 = dmn<2, int>;
  using SubDmn2 = dmn<4, double>;
  using SubDmn3 = dmn<3, unsigned int>;

  using ProductDmn = dca::func::dmn_variadic<dmn_0<SubDmn1>, dmn_0<SubDmn2>, dmn_0<SubDmn3>>;

  ProductDmn test_dmn;

  EXPECT_NO_THROW(test_dmn(0, 0, 0));
  EXPECT_NO_THROW(test_dmn(SubDmn1::get_size() - 1, 0, 0));
  EXPECT_NO_THROW(test_dmn(1, 3, 2));

  // Index too big.
  EXPECT_THROW(test_dmn(SubDmn1::get_size(), 0, 0), std::runtime_error);
  EXPECT_THROW(test_dmn(0, 4, 0), std::runtime_error);
  EXPECT_THROW(test_dmn(1, 0, 3), std::runtime_error);

  // Index too small.
  EXPECT_THROW(test_dmn(-1, 0, 0), std::runtime_error);
  EXPECT_THROW(test_dmn(1, -1, 0), std::runtime_error);
  EXPECT_THROW(test_dmn(0, 2, -3), std::runtime_error);
}

#endif  // NDEBUG
