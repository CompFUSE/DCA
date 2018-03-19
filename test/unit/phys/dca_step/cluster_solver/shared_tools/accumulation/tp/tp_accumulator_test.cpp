// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements a no-change test for the two point accumulation of a mock configuration.

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/tp_accumulator.hpp"

#include <array>
#include "gtest/gtest.h"
#include <string>
#include <vector>

#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/function/function_utils.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

struct ConfigElement {
  double get_tau() const {
    return tau_;
  }
  double get_left_band() const {
    return band_;
  }
  double get_right_band() const {
    return band_;
  }
  double get_left_site() const {
    return r_;
  }
  double get_right_site() const {
    return r_;
  }

  int band_;
  int r_;
  double tau_;
};
using Configuration = std::array<std::vector<ConfigElement>, 2>;
using MatrixPair = std::array<dca::linalg::Matrix<double, dca::linalg::CPU>, 2>;

using G0Setup = typename dca::testing::G0Setup<dca::testing::LatticeBilayer>;
const std::string input_dir =
    DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/";

void prepareRandomConfig(Configuration& config, MatrixPair& M, std::array<int, 2> n);

TEST_F(G0Setup, Accumulate) {
  const std::array<int, 2> n{18, 22};
  MatrixPair M;
  Configuration config;
  prepareRandomConfig(config, M, n);

  auto& parameters = G0Setup::parameters;
  auto& data = *G0Setup::data;
  using TpDomain = typename Data::TpGreensFunction::this_domain_type;
  dca::func::function<std::complex<double>, TpDomain> G4_check;

  dca::io::HDF5Reader reader;

  reader.open_file(input_dir + "G4_result.hdf5");

  for (const dca::phys::FourPointType type :
       {dca::phys::PARTICLE_HOLE_TRANSVERSE, dca::phys::PARTICLE_HOLE_MAGNETIC,
        dca::phys::PARTICLE_HOLE_CHARGE, dca::phys::PARTICLE_PARTICLE_UP_DOWN}) {
    parameters.set_four_point_type(type);
    dca::phys::solver::accumulator::TpAccumulator<Parameters> accumulator(
        data.G0_k_w_cluster_excluded, parameters);

    const int sign = 1;
    accumulator.accumulate(M, config, sign);

    G4_check.set_name("G4_" + toString(type));
    reader.execute(G4_check);
    const auto diff = dca::func::utils::difference(accumulator.get_sign_times_G4(), G4_check);
    EXPECT_GT(5e-7, diff.l_inf);
  }

  reader.close_file();
}

TEST_F(G0Setup, SumTo) {
  auto& parameters = G0Setup::parameters;
  auto& data = *G0Setup::data;
  parameters.set_four_point_type(dca::phys::PARTICLE_HOLE_TRANSVERSE);

  using Accumulator =
      dca::phys::solver::accumulator::TpAccumulator<G0Setup::Parameters, dca::linalg::CPU>;
  Accumulator accumulator1(data.G0_k_w_cluster_excluded, parameters);
  Accumulator accumulator2(data.G0_k_w_cluster_excluded, parameters);
  Accumulator accumulator3(data.G0_k_w_cluster_excluded, parameters);

  const std::array<int, 2> n{4, 4};
  const int sign = -1;
  MatrixPair M1, M2;
  Configuration config1, config2;
  prepareRandomConfig(config1, M1, n);
  prepareRandomConfig(config2, M2, n);

  accumulator1.accumulate(M1, config1, sign);
  accumulator2.accumulate(M2, config2, sign);

  accumulator3.accumulate(M1, config1, sign);
  accumulator3.accumulate(M2, config2, sign);

  accumulator1.sumTo(accumulator2);
  const auto diff = dca::func::utils::difference(accumulator3.get_sign_times_G4(),
                                                 accumulator2.get_sign_times_G4());
  EXPECT_GT(5e-7, diff.l_inf);
}

void prepareRandomConfig(Configuration& config, MatrixPair& M, const std::array<int, 2> ns) {
  dca::math::random::StdRandomWrapper<std::ranlux48_base> rng(0, 1, 0);

  for (int s = 0; s < 2; ++s) {
    const int n = ns[s];
    config[s].resize(n);
    M[s].resize(n);
    for (int i = 0; i < n; ++i) {
      const double tau = rng() - 0.5;
      const int r = rng() * G0Setup::RDmn::dmn_size();
      const int b = rng() * G0Setup::BDmn::dmn_size();
      config[s][i] = ConfigElement{b, r, tau};
    }

    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i)
        M[s](i, j) = 2 * rng() - 1.;
  }
}
