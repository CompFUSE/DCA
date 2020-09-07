// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Test for phase.hpp

#include "dca/math/util/phase.hpp"

#include "gtest/gtest.h"

template <typename T>
class PhaseTest : public ::testing::Test {};

using Types = ::testing::Types<std::complex<double>, int>;
TYPED_TEST_CASE(PhaseTest, Types);

TEST(PhaseTest, MultiplyCmplx) {
  using Complex = std::complex<double>;

  const Complex a{1, 3};
  const Complex b{4, -6};

  dca::math::Phase<Complex> phase1;
  phase1.multiply(a);
  phase1.multiply(b);

  dca::math::Phase<Complex> phase2(a * b);

  EXPECT_NEAR(std::abs(phase1.getSign() - phase2.getSign()), 0, 1e-14);

  phase2.multiply(phase2);
  EXPECT_NEAR(std::abs(std::pow(phase1.getSign(), 2) - phase2.getSign()), 0, 1e-14);
}

TEST(PhaseTest, DivideCmplx) {
  using Complex = std::complex<double>;

  const Complex a{0, -1};
  const Complex b{-3, 1.5};

  dca::math::Phase<Complex> phase1;
  phase1.multiply(a);
  phase1.divide(b);

  dca::math::Phase<Complex> phase2(a / b);

  EXPECT_NEAR(std::abs(phase1.getSign() - phase2.getSign()), 0, 1e-14);

  phase2.divide(phase2);
  phase1.reset();
  EXPECT_NEAR(std::abs(phase1.getSign() - phase2.getSign()), 0, 1e-14);
}

TYPED_TEST(PhaseTest, Null) {
  using Scalar = TypeParam;

  dca::math::Phase<Scalar> phase1;
  phase1.multiply(Scalar(0));

  dca::math::Phase<Scalar> phase2;
  phase2.makeNull();

  // Further operations are ignored until a reset.
  phase1.multiply(Scalar(1));
  phase2.divide(phase1);

  EXPECT_NEAR(std::abs(phase1.getSign() - Scalar(0)), 0, 1e-14);
  EXPECT_NEAR(std::abs(phase2.getSign() - Scalar(0)), 0, 1e-14);
}
