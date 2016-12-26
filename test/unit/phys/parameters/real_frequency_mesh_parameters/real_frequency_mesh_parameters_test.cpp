// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests real_frequency_mesh_parameters.hpp
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/phys/parameters/real_frequency_mesh_parameters.hpp"
#include "gtest/gtest.h"
#include "dca/io/json/json_reader.hpp"

TEST(RealFrequencyMeshParametersTest, DefaultValues) {
  dca::phys::params::RealFrequencyMeshParameters pars;

  EXPECT_EQ(-10., pars.get_min_real_frequency());
  EXPECT_EQ(10., pars.get_max_real_frequency());
  EXPECT_EQ(128, pars.get_real_frequencies());
  EXPECT_EQ(0.01, pars.get_imaginary_damping());
}

TEST(RealFrequencyMeshParametersTest, ReadAll) {
  dca::io::JSONReader reader;
  dca::phys::params::RealFrequencyMeshParameters pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/real_frequency_mesh_parameters/input_read_all.json");
  pars.readWrite(reader);
  reader.close_file();

  EXPECT_EQ(-8., pars.get_min_real_frequency());
  EXPECT_EQ(8., pars.get_max_real_frequency());
  EXPECT_EQ(64, pars.get_real_frequencies());
  EXPECT_EQ(0.001, pars.get_imaginary_damping());
}
