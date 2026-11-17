// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Andrei Plamada (plamada@phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests mpi_collective_sum.hpp.
//

#include "dca/parallel/mpi_concurrency/mpi_collective_sum.hpp"

#include <cmath>
#include <complex>
#include <numeric>
#include <vector>

#include "gtest/gtest.h"

#include "dca/testing/minimalist_printer.hpp"

class MPICollectiveSumTest : public ::testing::Test {
protected:
  MPICollectiveSumTest() {
    size_ = sum_interface_.get_size();
    rank_ = sum_interface_.get_id();
  }

  int size_;
  int rank_;
  constexpr static double epsilon_ = 1e-14;

  dca::parallel::MPICollectiveSum sum_interface_;
};

TEST_F(MPICollectiveSumTest, SumScalar) {
  int scalar_test;
  int scalar_expected;

  scalar_test = rank_;
  sum_interface_.sum(scalar_test);

  scalar_expected = size_ * (size_ - 1) / 2;

  MPI_Barrier(MPI_COMM_WORLD);

  EXPECT_EQ(scalar_expected, scalar_test);
}

TEST_F(MPICollectiveSumTest, SumFunction) {
  using TestDomain = dca::func::dmn_0<dca::func::dmn<2, int>>;

  dca::func::function<double, TestDomain> function_test("test");
  dca::func::function<double, TestDomain> function_expected("expected");

  for (int i = 0; i < function_test.size(); i++)
    function_test(i) = i * rank_;

  sum_interface_.sum(function_test);
  
  
  for (int i = 0; i < function_test.size(); i++)
    function_expected(i) = i * size_ * (size_ - 1) / 2;

  MPI_Barrier(MPI_COMM_WORLD);

  for (int i = 0; i < function_test.size(); i++)
    EXPECT_EQ(function_expected(i), function_test(i));
}

TEST_F(MPICollectiveSumTest, DelayedLeaveOneOutSumScalars) {
  const double real_local = 1.25 + 2. * rank_;
  const std::complex<double> complex_local(0.5 + rank_, 2. - 0.75 * rank_);
  double real_immediate = real_local;
  std::complex<double> complex_immediate = complex_local;
  double real_delayed = real_local;
  std::complex<double> complex_delayed = complex_local;
  double ordinary_delayed = rank_ + 1.;

  sum_interface_.leaveOneOutSum(real_immediate);
  sum_interface_.leaveOneOutSum(complex_immediate);
  sum_interface_.leaveOneOutSum(real_delayed, true);
  sum_interface_.leaveOneOutSum(complex_delayed, true);
  sum_interface_.delayedSum(ordinary_delayed);
  sum_interface_.resolveSums();

  double real_expected = size_ == 1 ? real_local : 0.;
  std::complex<double> complex_expected = size_ == 1 ? complex_local : std::complex<double>{};
  double ordinary_expected = 0.;
  for (int rank = 0; rank < size_; ++rank) {
    if (rank != rank_) {
      real_expected += 1.25 + 2. * rank;
      complex_expected += std::complex<double>(0.5 + rank, 2. - 0.75 * rank);
    }
    ordinary_expected += rank + 1.;
  }

  EXPECT_DOUBLE_EQ(real_expected, real_immediate);
  EXPECT_DOUBLE_EQ(complex_expected.real(), complex_immediate.real());
  EXPECT_DOUBLE_EQ(complex_expected.imag(), complex_immediate.imag());
  EXPECT_DOUBLE_EQ(real_expected, real_delayed);
  EXPECT_DOUBLE_EQ(complex_expected.real(), complex_delayed.real());
  EXPECT_DOUBLE_EQ(complex_expected.imag(), complex_delayed.imag());
  EXPECT_DOUBLE_EQ(ordinary_expected, ordinary_delayed);
}

TEST_F(MPICollectiveSumTest, DelayedLeaveOneOutSumFunctions) {
  using TestDomain = dca::func::dmn_0<dca::func::dmn<3, int>>;
  using RealFunction = dca::func::function<double, TestDomain>;
  using ComplexFunction = dca::func::function<std::complex<double>, TestDomain>;

  RealFunction real_immediate("real-immediate");
  RealFunction real_delayed("real-delayed");
  RealFunction ordinary_delayed("ordinary-delayed");
  ComplexFunction complex_immediate("complex-immediate");
  ComplexFunction complex_delayed("complex-delayed");
  for (int i = 0; i < real_immediate.size(); ++i) {
    real_immediate(i) = 1. + rank_ * (i + 1.);
    real_delayed(i) = real_immediate(i);
    ordinary_delayed(i) = real_immediate(i);
    complex_immediate(i) = {2. + rank_ * (i + 1.), -1. + rank_ * (2. * i + 1.)};
    complex_delayed(i) = complex_immediate(i);
  }

  sum_interface_.leaveOneOutSum(real_immediate);
  sum_interface_.leaveOneOutSum(complex_immediate);
  sum_interface_.leaveOneOutSum(real_delayed, true);
  sum_interface_.leaveOneOutSum(complex_delayed, true);
  sum_interface_.delayedSum(ordinary_delayed);
  sum_interface_.resolveSums();

  for (int i = 0; i < real_immediate.size(); ++i) {
    // The single-rank API intentionally leaves the local contribution unchanged.
    double real_expected = size_ == 1 ? 1. + rank_ * (i + 1.) : 0.;
    std::complex<double> complex_expected =
        size_ == 1 ? std::complex<double>(2. + rank_ * (i + 1.), -1. + rank_ * (2. * i + 1.))
                   : std::complex<double>{};
    for (int rank = 0; rank < size_; ++rank) {
      if (rank != rank_) {
        real_expected += 1. + rank * (i + 1.);
        complex_expected += std::complex<double>(2. + rank * (i + 1.),
                                                 -1. + rank * (2. * i + 1.));
      }
    }
    const double local = 1. + rank_ * (i + 1.);
    EXPECT_DOUBLE_EQ(size_ == 1 ? local : real_expected + local, ordinary_delayed(i));
    EXPECT_DOUBLE_EQ(real_expected, real_immediate(i));
    EXPECT_DOUBLE_EQ(real_expected, real_delayed(i));
    EXPECT_DOUBLE_EQ(complex_expected.real(), complex_immediate(i).real());
    EXPECT_DOUBLE_EQ(complex_expected.imag(), complex_immediate(i).imag());
    EXPECT_DOUBLE_EQ(complex_expected.real(), complex_delayed(i).real());
    EXPECT_DOUBLE_EQ(complex_expected.imag(), complex_delayed(i).imag());
  }
}

TEST_F(MPICollectiveSumTest, LeaveOneOutAvgAndSum) {
  std::vector<double> values(size_);
  double sum = 0.;
  for (int i = 0; i < size_; ++i) {
    values[i] = 3.14 + i;
    sum += values[i];
  }
  // These APIs intentionally leave the local value unchanged for one rank.
  const double expected = size_ == 1 ? values[rank_] : (sum - values[rank_]) / double(size_ - 1);
  const double expected_sum = size_ == 1 ? values[rank_] : expected * (size_ - 1);

  // Check scalar version.
  double scalar = values[rank_];
  double sum_one_out = scalar;
  sum_interface_.leaveOneOutAvg(scalar);
  sum_interface_.leaveOneOutSum(sum_one_out);

  MPI_Barrier(MPI_COMM_WORLD);

  EXPECT_DOUBLE_EQ(expected, scalar);
  EXPECT_DOUBLE_EQ(expected_sum, sum_one_out);

  // Check dca::func::function version.
  using TestDomain = dca::func::dmn_0<dca::func::dmn<2, int>>;
  dca::func::function<double, TestDomain> f;
  f(0) = values[rank_];
  f(1) = 0.;

  dca::func::function<double, TestDomain> f_sum(f);
  sum_interface_.leaveOneOutAvg(f);
  sum_interface_.leaveOneOutSum(f_sum);

  MPI_Barrier(MPI_COMM_WORLD);

  EXPECT_DOUBLE_EQ(expected, f(0));
  EXPECT_DOUBLE_EQ(expected_sum, f_sum(0));

  EXPECT_DOUBLE_EQ(0., f(1));
  EXPECT_DOUBLE_EQ(0., f_sum(1));
}

TEST_F(MPICollectiveSumTest, JackknifeErrorReal) {
  using TestDomain = dca::func::dmn_0<dca::func::dmn<3, int>>;
  using FunctionType = dca::func::function<double, TestDomain>;
  FunctionType original;
  if (size_ == 1) {
    original(0) = 7.;
    original(1) = 9.;
    original(2) = 11.;
    for (bool overwrite : {false, true}) {
      auto replicas = original;
      const auto error = sum_interface_.jackknifeError(replicas, overwrite);
      for (int i = 0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ(0., error(i));
        EXPECT_DOUBLE_EQ(original(i), replicas(i));
      }
    }
    return;
  }

  // Explicitly delete each rank's observation, then evaluate the estimator.
  // The ratio uses unequal denominators, so its replica average is not the pooled estimate.
  std::vector<double> means(size_), ratios(size_);
  double numerator = 0., denominator = 0.;
  for (int r = 0; r < size_; ++r) {
    means[r] = (size_ * (size_ + 1.) / 2. - (r + 1.)) / (size_ - 1.);
    numerator += (r + 1.) * (r + 1.);
    denominator += r + 1.;
  }
  for (int r = 0; r < size_; ++r)
    ratios[r] = (numerator - (r + 1.) * (r + 1.)) / (denominator - (r + 1.));

  // Independent oracle: sum_{r,s}(t_r-t_s)^2 = 2*N*sum_r(t_r-t_avg)^2.
  double pairwise = 0.;
  for (double a : ratios)
    for (double b : ratios)
      pairwise += (a - b) * (a - b);
  const double ratio_se = std::sqrt((size_ - 1.) * pairwise / (2. * size_ * size_));
  const double ratio_average = std::accumulate(ratios.begin(), ratios.end(), 0.) / size_;
  ASSERT_GT(std::abs(numerator / denominator - ratio_average), 1e-6);
  original(0) = means[rank_];
  original(1) = 9.;
  original(2) = ratios[rank_];

  for (bool overwrite : {false, true}) {
    auto replicas = original;
    const auto error = sum_interface_.jackknifeError(replicas, overwrite);
    // Analytic standard error of the mean of raw observations 1..N.
    EXPECT_NEAR(std::sqrt((size_ + 1.) / 12.), error(0), epsilon_);
    EXPECT_DOUBLE_EQ(0., error(1));
    EXPECT_NEAR(ratio_se, error(2), epsilon_);
    if (!overwrite) {
      for (int i = 0; i < 3; ++i)
        EXPECT_DOUBLE_EQ(original(i), replicas(i));
    }
    else {
      EXPECT_NEAR((size_ + 1.) / 2., replicas(0), epsilon_);
      EXPECT_DOUBLE_EQ(9., replicas(1));
      EXPECT_NEAR(ratio_average, replicas(2), epsilon_);
    }
  }
}

TEST_F(MPICollectiveSumTest, JackknifeErrorComplex) {
  using TestDomain = dca::func::dmn_0<dca::func::dmn<3, int>>;
  using FunctionType = dca::func::function<std::complex<double>, TestDomain>;
  FunctionType original;
  if (size_ == 1) {
    original(0) = {7., -2.};
    original(1) = {4., -3.};
    original(2) = {8., 5.};
    for (bool overwrite : {false, true}) {
      auto replicas = original;
      const auto error = sum_interface_.jackknifeError(replicas, overwrite);
      for (int i = 0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ(0., error(i).real());
        EXPECT_DOUBLE_EQ(0., error(i).imag());
        EXPECT_EQ(original(i), replicas(i));
      }
    }
    return;
  }

  // Different real/imaginary observations test marginal errors independently.
  // Also construct nonlinear ratio replicas rather than passing raw observations.
  std::vector<std::complex<double>> observations(size_), means(size_), ratios(size_);
  std::complex<double> total{}, numerator{}, denominator{};
  for (int r = 0; r < size_; ++r) {
    observations[r] = {r + 1., (r + 1.) * (r + 1.)};
    total += observations[r];
    const std::complex<double> x(r + 1., .5 * (r + 1.));
    numerator += x * x;
    denominator += x;
  }
  for (int r = 0; r < size_; ++r) {
    means[r] = (total - observations[r]) / double(size_ - 1);
    const std::complex<double> x(r + 1., .5 * (r + 1.));
    ratios[r] = (numerator - x * x) / (denominator - x);
  }
  const auto pairwise_error = [this](const auto& samples, bool imaginary) {
    double pairwise = 0.;
    for (const auto& a : samples)
      for (const auto& b : samples) {
        const double d = imaginary ? a.imag() - b.imag() : a.real() - b.real();
        pairwise += d * d;
      }
    return std::sqrt((size_ - 1.) * pairwise / (2. * size_ * size_));
  };
  std::complex<double> ratio_average{};
  for (const auto& ratio : ratios)
    ratio_average += ratio / double(size_);
  original(0) = means[rank_];
  original(1) = {4., -3.};
  original(2) = ratios[rank_];

  for (bool overwrite : {false, true}) {
    auto replicas = original;
    const auto error = sum_interface_.jackknifeError(replicas, overwrite);
    EXPECT_NEAR(std::sqrt((size_ + 1.) / 12.), error(0).real(), epsilon_);
    EXPECT_NEAR(pairwise_error(means, true), error(0).imag(), epsilon_);
    EXPECT_NE(error(0).real(), error(0).imag());
    EXPECT_DOUBLE_EQ(0., error(1).real());
    EXPECT_DOUBLE_EQ(0., error(1).imag());
    EXPECT_NEAR(pairwise_error(ratios, false), error(2).real(), epsilon_);
    EXPECT_NEAR(pairwise_error(ratios, true), error(2).imag(), epsilon_);
    if (!overwrite) {
      for (int i = 0; i < 3; ++i)
        EXPECT_EQ(original(i), replicas(i));
    }
    else {
      EXPECT_NEAR((size_ + 1.) / 2., replicas(0).real(), epsilon_);
      EXPECT_NEAR((size_ + 1.) * (2. * size_ + 1.) / 6., replicas(0).imag(), epsilon_);
      EXPECT_EQ(std::complex<double>(4., -3.), replicas(1));
      EXPECT_NEAR(ratio_average.real(), replicas(2).real(), epsilon_);
      EXPECT_NEAR(ratio_average.imag(), replicas(2).imag(), epsilon_);
    }
  }
}

TEST_F(MPICollectiveSumTest, ComputeCovarianceScalar) {
  using FunctionDomain = dca::func::dmn_0<dca::func::dmn<4, int>>;
  using CovarianceDomain = dca::func::dmn_variadic<FunctionDomain, FunctionDomain>;

  dca::func::function<double, FunctionDomain> f("f");
  dca::func::function<double, FunctionDomain> f_mean("f_mean");

  for (int i = 0; i < f.size(); ++i) {
    f(i) = i * rank_;
    f_mean(i) = f(i);
  }
  sum_interface_.sum(f_mean);
  f_mean /= size_;

  /// Calculate the equivalent covariance matrix without using function's etc
  std::vector<std::vector<double>> cov;
  for (int i = 0; i < f.size(); ++i)
    cov.emplace_back(f.size(), 0.0);

  auto covariance = [&cov](auto& samples, auto& f_estimated) {
    for (int i = 0; i < samples.size(); ++i)
      for (int j = 0; j < samples.size(); ++j)
        cov.at(i)[j] += (samples[i] - f_estimated[i]) * (samples[j] - f_estimated[j]);
  };

  std::vector<double> vec_means(f.size(), 0.0);

  for (int i = 0; i < f.size(); ++i)
    vec_means[i] = f_mean(i);

  std::vector<double> samples(f.size(), 0.0);
  for (int j = 0; j < size_; ++j) {
    for (int i = 0; i < f.size(); ++i) {
      samples[i] = i * j;
    }
    covariance(samples, vec_means);
  }
  for (int i = 0; i < samples.size(); ++i)
    for (int j = 0; j < samples.size(); ++j)
      cov.at(i)[j] /= size_;

  dca::func::function<double, CovarianceDomain> covariance_from_computeCovariance(
      "computeCovariance");
  dca::func::function<double, CovarianceDomain> covariance_from_computeCovarianceAndAvg(
      "computeCovarianceAndAvg");

  // Compute covariance matrix with respect to the *precomputed* mean of f.
  sum_interface_.computeCovariance(f, f_mean, covariance_from_computeCovariance);
  // Compute covariance matrix with respect to the mean of f.
  sum_interface_.computeCovarianceAndAvg(f, covariance_from_computeCovarianceAndAvg);

  MPI_Barrier(MPI_COMM_WORLD);

  for (int i = 0; i < f_mean.size(); ++i)
    EXPECT_DOUBLE_EQ(f_mean(i), f(i));

  for (int i = 0; i < f.size(); ++i) {
    for (int j = 0; j < f.size(); ++j) {
      EXPECT_DOUBLE_EQ(cov.at(j)[i], covariance_from_computeCovariance(i, j));
      EXPECT_DOUBLE_EQ(cov.at(j)[i], covariance_from_computeCovarianceAndAvg(i, j));
    }
  }
}

TEST_F(MPICollectiveSumTest, ComputeCovarianceComplex) {
  using FunctionDomain = dca::func::dmn_0<dca::func::dmn<2, int>>;
  using CovarianceDomain = dca::func::dmn_variadic<dca::func::dmn_0<dca::func::dmn<4, int>>,
                                                   dca::func::dmn_0<dca::func::dmn<4, int>>>;

  dca::func::function<std::complex<double>, FunctionDomain> f("f");
  dca::func::function<std::complex<double>, FunctionDomain> f_mean("f_mean");

  for (int i = 0; i < f.size(); ++i) {
    f(i) = std::complex<double>(rank_ * i, rank_ * (i + f.size()));
    f_mean(i) = f(i);
  }
  sum_interface_.sum(f_mean);

  MPI_Barrier(MPI_COMM_WORLD);

  f_mean /= size_;

  /// Calculate the equivalent covariance matrix without using function's etc
  /// The form of this matrix is quite bizarre that instead of being rank f complex it is instead
  /// rank 2f matrix containing the decomposed M_zz over M_xx, M_xy, M_yx, M_zz in its quadrants.

  /// Instead we will calculate the Covariance and Pseudo-covariance matrix which we use to check the M_xx, M_xy, etc.

  std::vector<std::vector<std::complex<double>>> cov;
  for (int i = 0; i < f.size(); ++i)
    cov.emplace_back(f.size(), std::complex<double>{0.0, 0.0});

  auto covariance = [&cov](auto& samples, auto& f_estimated) {
    for (int i = 0; i < samples.size(); ++i)
      for (int j = 0; j < samples.size(); ++j) {
        cov.at(i).at(j) += (samples[i] - f_estimated[i]) * std::conj(samples[j] - f_estimated[j]);
      }
  };

  std::vector<std::vector<std::complex<double>>> ps_cov;
  for (int i = 0; i < f.size(); ++i)
    ps_cov.emplace_back(f.size(), std::complex<double>{0.0, 0.0});

  auto ps_covariance = [&ps_cov](auto& samples, auto& f_estimated) {
    for (int i = 0; i < samples.size(); ++i)
      for (int j = 0; j < samples.size(); ++j)
        ps_cov.at(i).at(j) += (samples[i] - f_estimated[i]) * (samples[j] - f_estimated[j]);
  };

  std::vector<std::complex<double>> vec_means(f.size(), 0.0);

  for (int i = 0; i < f.size(); ++i)
    vec_means[i] = f_mean(i);

  std::vector<std::complex<double>> samples(f.size());
  for (int j = 0; j < size_; ++j) {
    for (int i = 0; i < f.size(); ++i) {
      samples[i] = std::complex<double>(j * i, j * (i + f.size()));
    }
    covariance(samples, vec_means);
    ps_covariance(samples, vec_means);
  }
  for (int i = 0; i < samples.size(); ++i) {
    for (int j = 0; j < samples.size(); ++j) {
      cov.at(i)[j] /= size_;
      ps_cov.at(i)[j] /= size_;
    }
  }
  dca::func::function<double, CovarianceDomain> covariance_from_computeCovariance(
      "computeCovariance");
  dca::func::function<double, CovarianceDomain> covariance_from_computeCovarianceAndAvg(
      "computeCovarianceAndAvg");

  // Compute covariance matrix with respect to the *precomputed* mean of f.
  sum_interface_.computeCovariance(f, f_mean, covariance_from_computeCovariance);
  MPI_Barrier(MPI_COMM_WORLD);

  // Compute covariance matrix with respect to the mean of f.
  sum_interface_.computeCovarianceAndAvg(f, covariance_from_computeCovarianceAndAvg);
  MPI_Barrier(MPI_COMM_WORLD);

  for (int i = 0; i < f_mean.size(); ++i) {
    EXPECT_DOUBLE_EQ(f_mean(i).real(), f(i).real());
    EXPECT_DOUBLE_EQ(f_mean(i).imag(), f(i).imag());
  }

  for (int i = 0; i < f.size(); ++i) {
    for (int j = 0; j < f.size(); ++j) {
      auto kzz_pzz_elem = std::complex<double>((cov.at(i)[j] + ps_cov.at(i)[j]));
      double real_real_elem = .5 * kzz_pzz_elem.real();
      kzz_pzz_elem = std::complex<double>((ps_cov.at(i)[j] - cov.at(i)[j]));
      double real_imag_elem = .5 * kzz_pzz_elem.imag();
      EXPECT_NEAR(real_real_elem, covariance_from_computeCovariance(i, j), 1.0E-12);
      EXPECT_NEAR(real_real_elem, covariance_from_computeCovarianceAndAvg(i, j), 1.0E-12);
      EXPECT_DOUBLE_EQ(real_imag_elem, covariance_from_computeCovariance(i, j + f.size()));
      EXPECT_DOUBLE_EQ(real_imag_elem, covariance_from_computeCovarianceAndAvg(i, j + f.size()));
    }
  }
}

TEST_F(MPICollectiveSumTest, AvgNormalizedMomenta) {
  using FunctionDomain = dca::func::dmn_0<dca::func::dmn<2, int>>;
  const std::vector<int> orders{3, 4};

  dca::func::function<double, FunctionDomain> f("f");

  std::vector<double> momenta(orders.size());

  for (int i = 0; i < f.size(); i++) {
    f(i) = rank_ * (i + 1);
  }

  momenta = sum_interface_.avgNormalizedMomenta(f, orders);
  MPI_Barrier(MPI_COMM_WORLD);

  if (size_ == 1) {
    // Normalized moments are undefined when the single-rank variance is zero.
    ASSERT_EQ(2, momenta.size());
    EXPECT_TRUE(std::isnan(momenta[0]));
    EXPECT_TRUE(std::isnan(momenta[1]));
    return;
  }

  std::vector<std::vector<double>> m_momenta;
  for (int i = 0; i < orders.size(); ++i)
    m_momenta.emplace_back(f.size(), 0.0);
  std::vector<double> samples(f.size(), 0.0);
  std::vector<double> var2(samples.size(), 0.0);
  auto unormalizedMomenta = [&var2, &m_momenta, &orders](auto& samples, auto& f_means) {
    for (int i = 0; i < samples.size(); i++) {
      const double diff = samples[i] - f_means[i];
      var2[i] += diff * diff;
      for (int iord = 0; iord < orders.size(); ++iord)
        m_momenta.at(iord)[i] += std::pow(diff, orders[iord]);
    }
  };

  std::vector<double> f_mean(f.size(), 0.0);
  for (int ir = 0; ir < size_; ++ir) {
    for (int i = 0; i < f.size(); ++i) {
      f_mean[i] += ir * (i + 1);
    }
  }

  for (int i = 0; i < f.size(); ++i)
    f_mean[i] /= size_;

  for (int ir = 0; ir < size_; ++ir) {
    for (int i = 0; i < f.size(); ++i) {
      samples[i] = ir * (i + 1);
    }

    unormalizedMomenta(samples, f_mean);
  }

  std::vector<double> expected_momenta_avg(orders.size(), 0.0);
  for (int i = 0; i < f.size(); ++i) {
    const double var = std::sqrt(var2[i] / size_);
    for (int j = 0; j < orders.size(); ++j)
      expected_momenta_avg[j] += std::abs(m_momenta.at(j)[i]) / (size_ * std::pow(var, orders[j]));
  }

  for (int ir = 0; ir < orders.size(); ++ir)
    expected_momenta_avg[ir] /= f.size();

  // Expected values obtained with python
  EXPECT_NEAR(0., momenta[0], 1e-8);
  EXPECT_NEAR(expected_momenta_avg[1], momenta[1], 1e-8);
}

TEST_F(MPICollectiveSumTest, DelayedSum) {
  using TestDomain1 = dca::func::dmn_0<dca::func::dmn<2, char>>;
  using TestDomain2 = dca::func::dmn_0<dca::func::dmn<3, char>>;

  dca::func::function<double, TestDomain1> function_test("test");
  dca::func::function<double, TestDomain1> function_expected("expected");

  dca::func::function<std::complex<double>, TestDomain2> function_cmplx_test("test");
  dca::func::function<std::complex<double>, TestDomain2> function_cmplx_expected("expected");

  for (int i = 0; i < function_test.size(); i++) {
    function_test(i) = i * rank_;
    function_cmplx_test(i) = std::complex<double>(i * rank_, 1);
    function_cmplx_expected(i) = std::complex<double>(i * size_ * (size_ - 1) / 2, size_);
    function_expected(i) = i * size_ * (size_ - 1) / 2;
  }

  double scalar_test = rank_;
  double scalar_expected = size_ * (size_ - 1) / 2;

  sum_interface_.delayedSum(function_test);
  sum_interface_.delayedSum(scalar_test);
  sum_interface_.delayedSum(function_cmplx_test);

  sum_interface_.resolveSums();
  MPI_Barrier(MPI_COMM_WORLD);

  EXPECT_EQ(scalar_expected, scalar_test);
  EXPECT_EQ(function_expected, function_test);
  EXPECT_EQ(function_cmplx_expected, function_cmplx_test);
}

int main(int argc, char** argv) {
  int result = 0;

  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  ::testing::InitGoogleTest(&argc, argv);

  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  if (rank != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();

  MPI_Finalize();

  return result;
}
