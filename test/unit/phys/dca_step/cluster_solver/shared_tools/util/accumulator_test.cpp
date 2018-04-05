// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests accumulator.hpp.

#include "dca/phys/dca_step/cluster_solver/shared_tools/util/accumulator.hpp"

#include <complex>
#include <stdexcept>

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
