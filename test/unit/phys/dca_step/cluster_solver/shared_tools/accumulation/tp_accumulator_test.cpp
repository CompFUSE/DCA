// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements a no-change test for the two point accumulation of a mock configuration.

#include "dca/phys/dca_step/cluster_solver/ctaux/accumulator/tp/accumulator_nonlocal_chi.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/accumulator/tp/accumulator_nonlocal_g.hpp"

#include <array>
#include <string>
#include <vector>
#include "gtest/gtest.h"

#include "dca/function/util/difference.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/four_point_type.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

constexpr bool update_baseline = false;

struct ConfigElement {
  double get_tau() const {
    return tau_;
  }
  double get_band() const {
    return band_;
  }
  double get_r_site() const {
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
    DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/";

void prepareRandomConfig(Configuration& config, MatrixPair& M, std::array<int, 2> n);

TEST_F(G0Setup, Accumulate) {
  const std::array<int, 2> n{18, 22};
  MatrixPair M;
  Configuration config;
  prepareRandomConfig(config, M, n);

  auto& parameters = G0Setup::parameters;
  auto& data = *G0Setup::data;

  Data::TpGreensFunction G4("G4");
  Data::ReducedTpGreensFunction G4_check(G4.get_name());

  dca::phys::solver::ctaux::accumulator_nonlocal_G<G0Setup::Parameters, G0Setup::Data> nonlocal_G_obj(
      parameters, data, 0);
  dca::phys::solver::ctaux::accumulator_nonlocal_chi<G0Setup::Parameters, G0Setup::Data> nonlocal_chi_obj(
      parameters, data, 0, G4);

  dca::io::HDF5Writer writer;
  dca::io::HDF5Reader reader;

  if (update_baseline)
    writer.open_file("tp_accumulator_test_baseline.hdf5");
  else
    reader.open_file(input_dir + "tp_accumulator_test_baseline.hdf5");

  for (const dca::phys::FourPointType type :
       {dca::phys::PARTICLE_HOLE_TRANSVERSE, dca::phys::PARTICLE_HOLE_MAGNETIC,
        dca::phys::PARTICLE_HOLE_CHARGE, dca::phys::PARTICLE_PARTICLE_UP_DOWN}) {
    G4 = 0;
    parameters.set_four_point_type(type);

    nonlocal_G_obj.execute(config[0], M[0], config[1], M[1]);
    const int sign = 1;
    nonlocal_chi_obj.execute(sign, nonlocal_G_obj);

    if (update_baseline) {
      G4.set_name("G4_" + toString(type));
      writer.execute(G4);
    }
    else {
      G4_check.set_name("G4_" + toString(type));
      reader.execute(G4_check);
      const auto diff = dca::func::util::difference(G4, G4_check);
      EXPECT_GT(1e-8, diff.l_inf);
    }
  }

  if (update_baseline)
    writer.close_file();
  else
    reader.close_file();
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
