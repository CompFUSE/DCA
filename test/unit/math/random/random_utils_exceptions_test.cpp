// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests the exceptions of the random number library utility functions.

#include "dca/math/random/random_utils.hpp"
#include <stdexcept>
#include "gtest/gtest.h"

using dca::math::random::detail::getGlobalId;
using dca::math::random::detail::generateSeed;

TEST(GetGlobalIdTest, InvalidLocalId) {
  EXPECT_THROW(getGlobalId(-1, 0, 1), std::logic_error);
}

TEST(GetGlobalIdTest, InvalidProcId) {
  EXPECT_THROW(getGlobalId(0, -1, 1), std::logic_error);
  EXPECT_THROW(getGlobalId(0, 1, 1), std::logic_error);
  EXPECT_THROW(getGlobalId(0, 1, 0), std::logic_error);
  EXPECT_THROW(getGlobalId(0, -1, 0), std::logic_error);
  EXPECT_THROW(getGlobalId(0, -7, -2), std::logic_error);
}

TEST(GenerateSeedTest, InvalidGlobalId) {
  EXPECT_THROW(generateSeed(-1), std::logic_error);
  EXPECT_THROW(generateSeed(-10, 42), std::logic_error);
}
