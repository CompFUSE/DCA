// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests cpe_parameters.hpp
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/phys/parameters/cpe_parameters.hpp"
#include "gtest/gtest.h"
#include "comp_library/IO_library/JSON/JSON.hpp"

TEST(CpeParametersTest, DefaultValues) {
  dca::phys::params::CpeParameters pars;

  EXPECT_FALSE(pars.do_CPE());
  EXPECT_EQ(64, pars.get_N_wn());
  EXPECT_EQ(1., pars.get_CPE_smoothing_factor());
  EXPECT_EQ(100, pars.get_max_CPE_iterations());
  EXPECT_EQ(0., pars.get_max_CPE_error());
  EXPECT_FALSE(pars.simulate_gaussian_noise());
  EXPECT_EQ(1, pars.get_nr_of_CPE_samples());
  EXPECT_EQ(0., pars.get_simulated_CPE_stddev());
  EXPECT_FALSE(pars.compute_free_spectrum());
  EXPECT_FALSE(pars.compute_lattice_spectrum());
  EXPECT_FALSE(pars.compute_cluster_spectrum());
}

TEST(CpeParametersTest, ReadAll) {
  IO::reader<IO::JSON> reader;
  dca::phys::params::CpeParameters pars;

  reader.open_file(DCA_SOURCE_DIR "/test/unit/phys/parameters/cpe_parameters/input_read_all.json");
  pars.readWrite(reader);
  reader.close_file();

  EXPECT_TRUE(pars.do_CPE());
  EXPECT_EQ(128, pars.get_N_wn());
  EXPECT_EQ(0.9, pars.get_CPE_smoothing_factor());
  EXPECT_EQ(200, pars.get_max_CPE_iterations());
  EXPECT_EQ(1.e-3, pars.get_max_CPE_error());
  EXPECT_TRUE(pars.simulate_gaussian_noise());
  EXPECT_EQ(10, pars.get_nr_of_CPE_samples());
  EXPECT_EQ(1.e-2, pars.get_simulated_CPE_stddev());
  EXPECT_TRUE(pars.compute_free_spectrum());
  EXPECT_TRUE(pars.compute_lattice_spectrum());
  EXPECT_TRUE(pars.compute_cluster_spectrum());
}
