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
// It checks the first random number drawn from a random number generator constructed with the
// default seed.

#include "dca/math/random/std_random_wrapper.hpp"
#include <vector>
#include "gtest/gtest.h"
#include "random_tests_helper.hpp"

namespace dca {
namespace testing {
// dca::testing::

// Hardcoded return value of the first call to the tested generators
template <typename Engine>
inline std::vector<double>& getFirstRandomNumbersDefaultSeed();
template <>
inline std::vector<double>& getFirstRandomNumbersDefaultSeed<std::mt19937>() {
  static std::vector<double> vec{0.42217690740692559, 0.96778114638316592, 0.23346981000143191,
                                 0.82172293527874019, 0.71580268548762127};
  return vec;
}
template <>
inline std::vector<double>& getFirstRandomNumbersDefaultSeed<std::mt19937_64>() {
  static std::vector<double> vec{0.58960821903135685, 0.027846668667417744, 0.96174772617893955,
                                 0.79004886762619775, 0.73670732232143421};
  return vec;
}
template <>
inline std::vector<double>& getFirstRandomNumbersDefaultSeed<std::ranlux48_base>() {
  static std::vector<double> vec{0.72777987290424495, 0.6717218115467265, 0.67916691101138849,
                                 0.040260746353792579, 0.94799739988482223};
  return vec;
}
template <>
inline std::vector<double>& getFirstRandomNumbersDefaultSeed<std::ranlux48>() {
  static std::vector<double> vec{0.72777987290424495, 0.6717218115467265, 0.67916691101138849,
                                 0.040260746353792579, 0.94799739988482223};
  return vec;
}

}  // testing
}  // dca

template <typename T>
class StdRandomWrapperTest : public ::testing::Test {};

using Engines = ::testing::Types<std::mt19937, std::mt19937_64, std::ranlux48_base, std::ranlux48>;
TYPED_TEST_CASE(StdRandomWrapperTest, Engines);

TYPED_TEST(StdRandomWrapperTest, DefaultSeed) {
  dca::math::random::StdRandomWrapper<TypeParam> rng(1, 4);

  EXPECT_EQ(1, rng.getGlobalId());
  EXPECT_EQ(0, rng.getInitialSeed());
  EXPECT_EQ(14520465407477959144u, rng.getSeed());

  const int checks = 5;
  for (int i = 0; i < checks; ++i)
    EXPECT_DOUBLE_EQ(dca::testing::getFirstRandomNumbersDefaultSeed<TypeParam>()[i], rng());

  dca::testing::inUnitInterval(rng, 1000);
}
