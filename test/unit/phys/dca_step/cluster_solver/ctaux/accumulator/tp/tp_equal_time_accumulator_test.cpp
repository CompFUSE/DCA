// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific
// publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests the summation of equal time two particle results.

#include "dca/phys/dca_step/cluster_solver/ctaux/accumulator/tp/tp_equal_time_accumulator.hpp"

#include "gtest/gtest.h"

#include "dca/function/util/difference.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/accumulation_test.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

constexpr char input_file[] =
    DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/ctaux/accumulator/tp/input.json";

using TpEqualTimeAccumulatorTest =
    dca::testing::G0Setup<dca::testing::LatticeBilayer, dca::phys::solver::CT_AUX, input_file>;

using Configuration = std::array<std::vector<dca::phys::solver::ctaux::vertex_singleton>, 2>;
using Sample = std::array<dca::linalg::Matrix<double, dca::linalg::CPU>, 2>;

void buildConfiguration(Configuration& config, Sample& sample, int b_size, int r_size, double beta,
                        std::array<int, 2> n);

TEST_F(TpEqualTimeAccumulatorTest, AccumulateAndSum) {
  const std::array<int, 2> n{7, 5};
  Sample M1, M2;
  Configuration config1, config2;

  buildConfiguration(config1, M1, TpEqualTimeAccumulatorTest::BDmn::dmn_size(),
                     TpEqualTimeAccumulatorTest::RDmn::dmn_size(), parameters_.get_beta(), n);
  buildConfiguration(config2, M2, TpEqualTimeAccumulatorTest::BDmn::dmn_size(),
                     TpEqualTimeAccumulatorTest::RDmn::dmn_size(), parameters_.get_beta(), n);

  using Accumulator =
      dca::phys::solver::ctaux::TpEqualTimeAccumulator<G0Setup::Parameters, G0Setup::Data, double>;

  Accumulator accumulator_sum(parameters_, *data_, 0);
  Accumulator accumulator1(parameters_, *data_, 1);
  Accumulator accumulator2(parameters_, *data_, 2);
  Accumulator accumulator3(parameters_, *data_, 3);
  const int sign = -1;

  accumulator1.accumulateAll(config1[0], M1[0], config1[1], M1[1], sign);
  accumulator2.accumulateAll(config2[0], M2[0], config2[1], M2[1], sign);
  accumulator1.sumTo(accumulator_sum);
  accumulator2.sumTo(accumulator_sum);
  accumulator_sum.finalize();

  accumulator3.accumulateAll(config1[0], M1[0], config1[1], M1[1], sign);
  accumulator3.accumulateAll(config2[0], M2[0], config2[1], M2[1], sign);
  accumulator3.finalize();

  auto expect_near = [](const auto& f1, const auto f2) {
    auto diff = dca::func::util::difference(f1, f2);
    EXPECT_GE(5e-7, diff.l_inf);

    auto zero(f1);
    zero = 0;
    auto norm = dca::func::util::difference(f1, zero);
    EXPECT_GE(norm.l_inf, 0);
  };

  expect_near(accumulator_sum.get_G_r_t(), accumulator3.get_G_r_t());
  expect_near(accumulator_sum.get_charge_cluster_moment(), accumulator3.get_charge_cluster_moment());
  expect_near(accumulator_sum.get_magnetic_cluster_moment(),
              accumulator3.get_magnetic_cluster_moment());
  expect_near(accumulator_sum.get_dwave_pp_correlator(), accumulator3.get_dwave_pp_correlator());
}

void buildConfiguration(Configuration& config, Sample& sample, int b_size, int r_size, double beta,
                        std::array<int, 2> n) {
  static dca::math::random::StdRandomWrapper<std::ranlux48_base> rng(0, 1, 0);
  using namespace dca::phys;
  using namespace dca::phys::solver::ctaux;

  for (int s = 0; s < 2; ++s) {
    sample[s].resizeNoCopy(n[s]);

    const int band = rng() * b_size;
    const int r = rng() * r_size;
    const int dr = rng() * r_size;
    const double tau = rng() * beta;
    const auto field = rng() > 0.5 ? 1 : -1;

    for (int j = 0; j < n[s]; ++j) {
      auto vertex = dca::phys::solver::ctaux::vertex_singleton(
          band, s == 0 ? e_UP : e_DN, /*spin_orbitals.first*/ 0,
          /*spin_orbitals.second,=*/0, r, dr, tau, field ? HS_UP : HS_DN,
          field ? HS_FIELD_UP : HS_FIELD_DN, j);

      config[s].push_back(vertex);
      for (int i = 0; i < n[s]; ++i)
        sample[s](i, j) = 2 * rng() - 1;
    }
  }
}
