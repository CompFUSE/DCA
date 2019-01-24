// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file tests the reading and writing of the CT-AUX configuration to a buffer.

#include <dca/phys/dca_step/cluster_solver/cluster_solver_name.hpp>
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/read_write_config.hpp"

#include "gtest/gtest.h"

#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/cv.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

constexpr char input_name[] =
    DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/ctaux/structs/input.json";

using ReadWriteConfigTest =
    dca::testing::G0Setup<dca::testing::LatticeBilayer, dca::phys::solver::CT_AUX, input_name>;

TEST_F(ReadWriteConfigTest, All) {
  std::vector<double> random(100);
  for (auto& x : random)
    x = static_cast<double>(std::rand()) / RAND_MAX;
  Parameters::random_number_generator stub_rng(random);

  dca::phys::solver::ctaux::CV<Parameters>::get_H_interaction() = data_->H_interactions;

  dca::phys::solver::ctaux::CT_AUX_HS_configuration<Parameters> config(parameters_, stub_rng);
  config.initialize();

  dca::io::Buffer buffer;
  buffer << config;

  dca::phys::solver::ctaux::CT_AUX_HS_configuration<Parameters> config2(parameters_, stub_rng);
  buffer >> config2;

  EXPECT_EQ(config, config2);
}
