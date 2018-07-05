// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests the 2d 4-band model specialization of the ModelParameters class.
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/phys/parameters/model_parameters.hpp"
#include "gtest/gtest.h"
#include "dca/io/json/json_reader.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"

using PointGroup = dca::phys::domains::D4;

TEST(ModelParameters2d2bandTest, DefaultValues) {
  dca::phys::params::ModelParameters<
      dca::phys::models::TightBindingModel<dca::phys::models::fourband_lattice<PointGroup>>>
      pars;

  EXPECT_EQ(0., pars.get_ei0());
  EXPECT_EQ(0., pars.get_eb0());
  EXPECT_EQ(0., pars.get_t0());
  EXPECT_EQ(0., pars.get_ei1());
  EXPECT_EQ(0., pars.get_eb1());
  EXPECT_EQ(0., pars.get_t1());
  EXPECT_EQ(0., pars.get_U0());
  EXPECT_EQ(0., pars.get_U1());
  EXPECT_EQ(0., pars.get_V());
  EXPECT_EQ(0., pars.get_V_prime());
}

TEST(ModelParameters2d2bandTest, ReadAll) {
  dca::io::JSONReader reader;
  dca::phys::params::ModelParameters<
      dca::phys::models::TightBindingModel<dca::phys::models::fourband_lattice<PointGroup>>>
      pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/model_parameters/input_read_all_2d_4band.json");
  pars.readWrite(reader);
  reader.close_file();

  EXPECT_EQ(1., pars.get_ei0());
  EXPECT_EQ(2., pars.get_eb0());
  EXPECT_EQ(3., pars.get_t0());
  EXPECT_EQ(4., pars.get_ei1());
  EXPECT_EQ(5., pars.get_eb1());
  EXPECT_EQ(6., pars.get_t1());
  EXPECT_EQ(7., pars.get_U0());
  EXPECT_EQ(8., pars.get_U1());
  EXPECT_EQ(9., pars.get_V());
  EXPECT_EQ(10., pars.get_V_prime());
}
