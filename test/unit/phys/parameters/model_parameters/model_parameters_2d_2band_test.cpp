// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests the 2d 2-band model specialization of the ModelParameters class.
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/phys/parameters/model_parameters.hpp"
#include "gtest/gtest.h"
#include "comp_library/IO_library/JSON/JSON.hpp"
#include "phys_library/domains/cluster/symmetries/point_groups/2D/2D_square.h"

using PointGroup = D4;

TEST(ModelParameters2d2bandTest, DefaultValues) {
  dca::phys::params::ModelParameters<dca::phys::models::TightBindingModel<twoband_lattice<PointGroup>>> pars;

  EXPECT_EQ(0., pars.get_ei0());
  EXPECT_EQ(0., pars.get_eb0());
  EXPECT_EQ(0., pars.get_t0());
  EXPECT_EQ(0., pars.get_U0());
}

TEST(ModelParameters2d2bandTest, ReadAll) {
  IO::reader<IO::JSON> reader;
  dca::phys::params::ModelParameters<dca::phys::models::TightBindingModel<twoband_lattice<PointGroup>>> pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/model_parameters/input_read_all_2d_2band.json");
  pars.readWrite(reader);
  reader.close_file();

  EXPECT_EQ(1., pars.get_ei0());
  EXPECT_EQ(2., pars.get_eb0());
  EXPECT_EQ(3., pars.get_t0());
  EXPECT_EQ(4., pars.get_U0());
}
