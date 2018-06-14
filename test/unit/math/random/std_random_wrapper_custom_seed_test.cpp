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
// It checks the first random number drawn from a random number generator constructed with a custom
// seed.

#include "dca/math/random/std_random_wrapper.hpp"
#include "gtest/gtest.h"
#include "random_tests_helper.hpp"

namespace dca {
namespace testing {
// dca::testing::

// Hardcoded return value of the first call to the tested generators
template <typename Engine>
inline std::vector<double>& getFirstRandomNumbersCustomSeed();
template <>
inline std::vector<double>& getFirstRandomNumbersCustomSeed<std::mt19937>() {
  static std::vector<double> vec{0.82432286004943012, 0.64010055391130516, 0.25076816535662333,
                                 0.085330161194926341, 0.25730370259651308};
  return vec;
}
template <>
inline std::vector<double>& getFirstRandomNumbersCustomSeed<std::mt19937_64>() {
  static std::vector<double> vec{0.27940464731427572, 0.055566565448361416, 0.5621891456799587,
                                 0.41722249930558986, 0.095007484445949017};
  return vec;
}
template <>
inline std::vector<double>& getFirstRandomNumbersCustomSeed<std::ranlux48_base>() {
  static std::vector<double> vec{0.0068332189419341447, 0.86802894579254541, 0.25453437079295799,
                                 0.24682368736319102, 0.88335484565944544};
  return vec;
}
template <>
inline std::vector<double>& getFirstRandomNumbersCustomSeed<std::ranlux48>() {
  static std::vector<double> vec{0.0068332189419341447, 0.86802894579254541, 0.25453437079295799,
                                 0.24682368736319102, 0.88335484565944544};
  return vec;
}

}  // testing
}  // dca

template <typename T>
class StdRandomWrapperTest : public ::testing::Test {};

using Engines = ::testing::Types<std::mt19937, std::mt19937_64, std::ranlux48_base, std::ranlux48>;
TYPED_TEST_CASE(StdRandomWrapperTest, Engines);

TYPED_TEST(StdRandomWrapperTest, CustomSeed) {
  dca::math::random::StdRandomWrapper<TypeParam> rng(2, 4, 139);

  EXPECT_EQ(2, rng.getGlobalId());
  EXPECT_EQ(139, rng.getInitialSeed());
  EXPECT_EQ(15128777329869828874u, rng.getSeed());

  const int checks = 5;
  for (int i = 0; i < checks; ++i)
    EXPECT_DOUBLE_EQ(dca::testing::getFirstRandomNumbersCustomSeed<TypeParam>()[i], rng());

  dca::testing::inUnitInterval(rng, 1000);
}
