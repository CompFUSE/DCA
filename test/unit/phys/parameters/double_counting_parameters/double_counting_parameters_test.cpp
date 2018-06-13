// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests double_counting_parameters.hpp
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/phys/parameters/double_counting_parameters.hpp"
#include "gtest/gtest.h"
#include "dca/io/json/json_reader.hpp"

TEST(DoubleCountingParametersTest, DefaultValues) {
  dca::phys::params::DoubleCountingParameters pars;

  EXPECT_EQ("none", pars.get_double_counting_method());
  EXPECT_EQ(0., pars.get_double_counting_correction());
}

TEST(DoubleCountingParametersTest, ReadAll) {
  dca::io::JSONReader reader;
  dca::phys::params::DoubleCountingParameters pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/double_counting_parameters/input_read_all.json");
  pars.readWrite(reader);
  reader.close_file();

  EXPECT_EQ("constant-correction-with-U-correction", pars.get_double_counting_method());
  EXPECT_EQ(4.2, pars.get_double_counting_correction());
}

TEST(DoubleCountingParametersTest, ReadIllegal) {
  dca::io::JSONReader reader;
  dca::phys::params::DoubleCountingParameters pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/double_counting_parameters/input_read_illegal.json");
  EXPECT_THROW(pars.readWrite(reader), std::logic_error);
  reader.close_file();
}
