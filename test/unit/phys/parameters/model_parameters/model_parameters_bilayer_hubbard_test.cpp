// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests the square lattice bilayer Hubbard model specialization of the ModelParameters
// class.
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/platform/dca_gpu.h"
#include "dca/phys/parameters/model_parameters.hpp"
#include "gtest/gtest.h"
#include "dca/io/json/json_reader.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"

using PointGroup = dca::phys::domains::D4;

TEST(ModelParametersBilayerHubbardTest, DefaultValues) {
  dca::phys::params::ModelParameters<
      dca::phys::models::TightBindingModel<dca::phys::models::bilayer_lattice<PointGroup>>>
      pars;

  EXPECT_EQ(0., pars.get_t());
  EXPECT_EQ(0., pars.get_t_prime());
  EXPECT_EQ(0., pars.get_t_perp());
  EXPECT_EQ(0., pars.get_U());
  EXPECT_EQ(0., pars.get_V());
  EXPECT_EQ(0., pars.get_V_prime());
}

TEST(ModelParametersBilayerHubbardTest, ReadAll) {
  dca::io::JSONReader reader;
  dca::phys::params::ModelParameters<
      dca::phys::models::TightBindingModel<dca::phys::models::bilayer_lattice<PointGroup>>>
      pars;

  reader.open_file(
      DCA_SOURCE_DIR
      "/test/unit/phys/parameters/model_parameters/input_read_all_bilayer_hubbard.json");
  pars.readWrite(reader);
  reader.close_file();

  EXPECT_EQ(1., pars.get_t());
  EXPECT_EQ(0.5, pars.get_t_prime());
  EXPECT_EQ(0.2, pars.get_t_perp());
  EXPECT_EQ(8., pars.get_U());
  EXPECT_EQ(2., pars.get_V());
  EXPECT_EQ(2., pars.get_V_prime());
}
