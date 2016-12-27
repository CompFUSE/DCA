// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests four_point_parameters.hpp
//
// TODO: Add tests for get_buffer_size, pack, unpack, writing, and the computation of the 'exact'
//       momentum transfer and its index.

#include "dca/phys/parameters/four_point_parameters.hpp"
#include "gtest/gtest.h"
#include "dca/io/json/json_reader.hpp"

TEST(FourPointParametersTest, DefaultValues) {
  dca::phys::params::FourPointParameters<2> pars;

  std::vector<double> momentum_transfer_input_check{0., 0.};

  EXPECT_EQ(dca::phys::NONE, pars.get_four_point_type());
  EXPECT_EQ(momentum_transfer_input_check, pars.get_four_point_momentum_transfer_input());
  EXPECT_EQ(0, pars.get_four_point_frequency_transfer());
}

TEST(FourPointParametersTest, ReadAll) {
  dca::io::JSONReader reader;
  dca::phys::params::FourPointParameters<2> pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/four_point_parameters/input_read_all.json");
  pars.readWrite(reader);
  reader.close_file();

  std::vector<double> momentum_transfer_input_check{3.14, -1.57};

  EXPECT_EQ(dca::phys::PARTICLE_PARTICLE_SUPERCONDUCTING, pars.get_four_point_type());
  EXPECT_EQ(momentum_transfer_input_check, pars.get_four_point_momentum_transfer_input());
  EXPECT_EQ(1, pars.get_four_point_frequency_transfer());
}
