// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests the wrapper class for the random number library of C++11.
// It checks whether different random number generators get different seeds.

#include "dca/math/random/std_random_wrapper.hpp"
#include "gtest/gtest.h"

template <typename T>
class StdRandomWrapperTest : public ::testing::Test {};

using Engines = ::testing::Types<std::mt19937, std::mt19937_64, std::ranlux48_base, std::ranlux48>;
TYPED_TEST_CASE(StdRandomWrapperTest, Engines);

TYPED_TEST(StdRandomWrapperTest, UniqueSeedsWithDefaultSeed) {
  dca::math::random::StdRandomWrapper<TypeParam> rng_1(0, 2);
  dca::math::random::StdRandomWrapper<TypeParam> rng_2(0, 2);
  dca::math::random::StdRandomWrapper<TypeParam> rng_3(1, 2);

  EXPECT_NE(rng_1.getSeed(), rng_2.getSeed());
  EXPECT_NE(rng_1.getSeed(), rng_3.getSeed());
  EXPECT_NE(rng_2.getSeed(), rng_3.getSeed());
}

TYPED_TEST(StdRandomWrapperTest, UniqueSeedsWithCustomSeed) {
  dca::math::random::StdRandomWrapper<TypeParam> rng_4(0, 2, 77);
  dca::math::random::StdRandomWrapper<TypeParam> rng_5(0, 2, 77);
  dca::math::random::StdRandomWrapper<TypeParam> rng_6(1, 2, 77);

  EXPECT_NE(rng_4.getSeed(), rng_5.getSeed());
  EXPECT_NE(rng_4.getSeed(), rng_6.getSeed());
  EXPECT_NE(rng_5.getSeed(), rng_6.getSeed());
}
