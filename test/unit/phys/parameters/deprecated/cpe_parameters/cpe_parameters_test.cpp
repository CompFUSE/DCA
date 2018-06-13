// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests cpe_parameters.hpp
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/phys/parameters/cpe_parameters.hpp"
#include "gtest/gtest.h"
#include "dca/io/json/json_reader.hpp"

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
  dca::io::JSONReader reader;
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
