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
// This file tests the Rashba-Hubbard model specialization of the ModelParameters class.
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/platform/dca_gpu.h"
#include "dca/phys/parameters/model_parameters.hpp"
#include "dca/testing/gtest_h_w_warning_blocking.h"
#include "dca/io/json/json_reader.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/rashba_hubbard.hpp"

using Lattice = dca::phys::models::RashbaHubbard<dca::phys::domains::no_symmetry<2>>;

TEST(ModelParametersRanshbaHubbardTest, DefaultValues) {
  dca::phys::params::ModelParameters<dca::phys::models::TightBindingModel<Lattice>> pars;
  EXPECT_EQ(1., pars.get_t());
  EXPECT_EQ(0., pars.get_h());
  EXPECT_EQ(0., pars.get_lambda());
  EXPECT_EQ(0., pars.get_U());
}

TEST(ModelParametersSingleBandHubbardTest, ReadAll) {
  dca::io::JSONReader reader;
  dca::phys::params::ModelParameters<dca::phys::models::TightBindingModel<Lattice>> pars;

  reader.open_file(
      DCA_SOURCE_DIR
      "/test/unit/phys/parameters/model_parameters/input_read_all_rashba_hubbard.json");
  pars.readWrite(reader);
  reader.close_file();

  EXPECT_EQ(1., pars.get_t());
  EXPECT_EQ(0.5, pars.get_h());
  EXPECT_EQ(0.2, pars.get_lambda());
  EXPECT_EQ(8.0, pars.get_U());
}
