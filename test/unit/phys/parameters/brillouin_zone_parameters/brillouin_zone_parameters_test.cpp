// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests brillouin_zone_parameters.hpp
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/phys/parameters/brillouin_zone_parameters.hpp"
#include "gtest/gtest.h"
#include "comp_library/IO_library/JSON/JSON.hpp"

TEST(BrillouinZoneParametersTest, DefaultValues) {
  dca::phys::params::BrillouinZoneParameters pars;

  EXPECT_EQ("absolute", pars.get_coordinate_type());
  EXPECT_EQ(std::vector<std::string>(0), pars.get_coordinate_names());
  EXPECT_EQ(std::vector<std::vector<double>>(0), pars.get_Brillouin_zone_vectors());
}
