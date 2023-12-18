// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Peter W. Doak (doakpw@ornl.gov)
//
// This file tests accumulator.hpp.

#include "dca/phys/dca_step/cluster_solver/shared_tools/util/accumulator.hpp"
#include "dca/math/util/phase.hpp"

#include <complex>
#include <stdexcept>
#include <cmath>

#include "gtest/gtest.h"

TEST(AccumulatorTest, FloatingPointType) {
  dca::phys::solver::util::Accumulator<double> acc;

  EXPECT_EQ(0, acc.count());
  EXPECT_EQ(0., acc.sum());
  EXPECT_THROW(acc.mean(), std::logic_error);

  const int num_samples = 10;
  for (int i = 0; i < num_samples; ++i)
    acc.addSample(i);

  EXPECT_EQ(num_samples, acc.count());
  EXPECT_EQ(45, acc.sum());
  EXPECT_EQ(4.5, acc.mean());

  acc.reset();

  EXPECT_EQ(0, acc.count());
  EXPECT_EQ(0., acc.sum());
  EXPECT_THROW(acc.mean(), std::logic_error);

  const double d = 3.14;
  acc.addSample(d);

  EXPECT_EQ(1, acc.count());
  EXPECT_EQ(d, acc.sum());
  EXPECT_EQ(d, acc.mean());
}

TEST(AccumulatorTest, IntegerType) {
  dca::phys::solver::util::Accumulator<int> acc;

  EXPECT_EQ(0, acc.count());
  EXPECT_EQ(0, acc.sum());
  EXPECT_THROW(acc.mean(), std::logic_error);

  const int num_samples = 10;
  for (int i = 0; i < num_samples; ++i)
    acc.addSample(i);

  EXPECT_EQ(num_samples, acc.count());
  EXPECT_EQ(45, acc.sum());
  EXPECT_EQ(4.5, acc.mean());

  acc.reset();

  EXPECT_EQ(0, acc.count());
  EXPECT_EQ(0, acc.sum());
  EXPECT_THROW(acc.mean(), std::logic_error);

  const int j = 42;
  acc.addSample(j);

  EXPECT_EQ(1, acc.count());
  EXPECT_EQ(j, acc.sum());
  EXPECT_EQ(j, acc.mean());
}

TEST(AccumulatorTest, ComplexType) {
  dca::phys::solver::util::Accumulator<std::complex<float>> acc;

  EXPECT_EQ(0, acc.count());
  EXPECT_EQ(std::complex<float>(0.), acc.sum());
  EXPECT_THROW(acc.mean(), std::logic_error);

  const int num_samples = 10;
  for (int i = 0; i < num_samples; ++i)
    acc.addSample(std::complex<float>(i, i + 10));

  EXPECT_EQ(num_samples, acc.count());
  EXPECT_EQ(std::complex<float>(45, 145), acc.sum());
  EXPECT_EQ(std::complex<float>(4.5, 14.5), acc.mean());

  acc.reset();

  EXPECT_EQ(0, acc.count());
  EXPECT_EQ(std::complex<float>(0.), acc.sum());
  EXPECT_THROW(acc.mean(), std::logic_error);

  const std::complex<float> c(3.14, 2.72);
  acc.addSample(c);

  EXPECT_EQ(1, acc.count());
  EXPECT_EQ(c, acc.sum());
  EXPECT_EQ(c, acc.mean());
}

TEST(AccumulatorTest, Phase) {
  dca::phys::solver::util::Accumulator<dca::math::Phase<std::complex<float>>> acc;
  EXPECT_EQ(0, acc.count());
  dca::math::Phase<std::complex<float>> phase;

  constexpr double mag{1.0};
  phase.multiply(std::polar(mag, M_PI * 0.50));
  acc.addSample(phase.getSign());
  EXPECT_NEAR(1.0, std::imag(acc.sum()), 1E-4);

  dca::math::Phase<std::complex<float>> phase2;
  phase2.multiply(std::polar(mag, M_PI * 1));
  acc.addSample(phase2.getSign());
  EXPECT_NEAR(1.0, std::imag(acc.sum()), 1E-4);

  phase2.multiply(std::polar(mag, M_PI * 0.5));
  acc.addSample(phase2.getSign());
  EXPECT_NEAR(0.0, std::imag(acc.sum()), 1E-4);

  EXPECT_EQ(acc.count(), 3);
}
