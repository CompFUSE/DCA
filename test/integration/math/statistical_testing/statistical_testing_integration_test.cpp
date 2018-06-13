// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Integration test for the StatiscalTesting class. It tests the error detection power on mock data.

#include "dca/math/statistical_testing/statistical_testing.hpp"

#include "gtest/gtest.h"
#include <random>
#include <vector>

#include "dca/function/function.hpp"
#include "dca/function/domains/dmn.hpp"
#include "dca/function/domains/dmn_0.hpp"
#include "dca/function/domains/dmn_variadic.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/blas/blas2.hpp"
#include "dca/linalg/matrixop.hpp"

// Typical number of processes (independent samples) and independent degrees of freedom used in a
// validation test. See ctaux_square_lattice_validation.cpp.
constexpr int n_meas = 300;
constexpr int sample_dimension = 30;

using SampleDmn = dca::func::dmn_0<dca::func::dmn<sample_dimension>>;
using NMeasDmn = dca::func::dmn_0<dca::func::dmn<n_meas>>;
using Function = dca::func::function<double, SampleDmn>;
using Covariance = dca::func::function<double, dca::func::dmn_variadic<SampleDmn, SampleDmn>>;
using FunctionSamples = dca::func::function<double, dca::func::dmn_variadic<SampleDmn, NMeasDmn>>;

Covariance computeCovariance(const FunctionSamples& samples, const Function& avg);
FunctionSamples generateSamples();
Function average(const FunctionSamples& samples);

TEST(StatisticalTestingTest, Integration) {
  auto samples = generateSamples();
  Function avg = average(samples);
  Function expected;
  Covariance cov;

  auto computePvalue = [&](const double meas_err) {
    expected = meas_err;
    cov = computeCovariance(samples, expected);
    dca::math::StatisticalTesting test(avg, expected, cov);
    return test.computePValue(false, NMeasDmn::dmn_size(), true);
  };

  // The test pass if the measurements have no error.
  EXPECT_LT(0.05, computePvalue(0));

  // The test fails if the measurements have a systematic error.
  EXPECT_GE(0.05, computePvalue(1e-2));

  // The smallest error that can be detected with a 5% test is around 1e-6.
  EXPECT_LT(0.05, computePvalue(1e-7));
  EXPECT_GT(0.05, computePvalue(1e-6));
}

FunctionSamples generateSamples() {
  FunctionSamples samples;
  std::mt19937_64 rng(0);
  std::normal_distribution<double> gauss_distro(0., 1.);
  std::uniform_real_distribution<double> uniform_distro(-0.5, 0.5);

  // Generate a covariance matrix with variances of order of magnitude 1.
  dca::linalg::Matrix<double, dca::linalg::CPU> cov_sqrt(SampleDmn::dmn_size());
  dca::linalg::Matrix<double, dca::linalg::CPU> true_cov(SampleDmn::dmn_size());
  for (int j = 0; j < SampleDmn::dmn_size(); ++j)
    for (int i = 0; i < SampleDmn::dmn_size(); ++i)
      cov_sqrt(i, j) = uniform_distro(rng) + (i == j) * 1.;
  // Make sure that the covariance is positive definite.
  dca::linalg::matrixop::gemm('N', 'T', cov_sqrt, cov_sqrt, true_cov);
  //  true_cov.print();

  std::vector<double> sample(SampleDmn::dmn_size());
  for (int sample_id = 0; sample_id < SampleDmn::dmn_size(); ++sample_id) {
    for (int i = 0; i < SampleDmn::dmn_size(); ++i)
      sample[i] = gauss_distro(rng);
    // samples[i] <- true_cov * sample
    dca::linalg::blas::gemv("N", SampleDmn::dmn_size(), SampleDmn::dmn_size(), 1., true_cov.ptr(),
                            true_cov.leadingDimension(), sample.data(), 1, 0.,
                            &samples(0, sample_id), 1);
  }

  return samples;
}

Covariance computeCovariance(const FunctionSamples& samples, const Function& avg) {
  Covariance cov;

  for (int sample_id = 0; sample_id < SampleDmn::dmn_size(); ++sample_id)
    for (int j = 0; j < SampleDmn::dmn_size(); ++j)
      for (int i = 0; i < SampleDmn::dmn_size(); ++i)
        cov(i, j) += (samples(i, sample_id) - avg(i)) * (samples(j, sample_id) - avg(j));

  cov /= NMeasDmn::dmn_size();
  return cov;
}

Function average(const FunctionSamples& samples) {
  Function avg;
  for (int sample_id = 0; sample_id < SampleDmn::dmn_size(); ++sample_id)
    for (int i = 0; i < SampleDmn::dmn_size(); ++i)
      avg(i) += samples(i, sample_id);

  avg /= NMeasDmn::dmn_size();
  return avg;
}
