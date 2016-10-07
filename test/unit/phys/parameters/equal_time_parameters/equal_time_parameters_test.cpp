// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests equal_time_parameters.hpp
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/phys/parameters/equal_time_parameters.hpp"
#include "gtest/gtest.h"
#include "dca/io/json/json_reader.hpp"

TEST(EqualTimeParametersTest, DefaultValues) {
  dca::phys::params::EqualTimeParameters pars;
  EXPECT_FALSE(pars.do_equal_time_measurements());
}

TEST(EqualTimeParametersTest, ReadAll) {
  dca::io::JSONReader reader;
  dca::phys::params::EqualTimeParameters pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/equal_time_parameters/input_read_all.json");
  pars.readWrite(reader);
  reader.close_file();

  EXPECT_TRUE(pars.do_equal_time_measurements());
}
