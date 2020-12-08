// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests MCI_parameters.hpp
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/phys/parameters/mci_parameters.hpp"

#include <algorithm>  // for std::sort and std::unique
#include <fstream>
#include <limits>
#include <vector>

#include "gtest/gtest.h"

#include "dca/io/json/json_reader.hpp"

TEST(MciParametersTest, DefaultValues) {
  dca::phys::params::MciParameters pars;

  EXPECT_EQ(985456376, pars.get_seed());
  EXPECT_EQ(20, pars.get_warm_up_sweeps());
  EXPECT_EQ(std::vector<double>{1.}, pars.get_sweeps_per_measurement());
  EXPECT_EQ(std::vector<int>{100}, pars.get_measurements());
  EXPECT_EQ(dca::phys::ErrorComputationType::NONE, pars.get_error_computation_type());
  EXPECT_EQ(1, pars.get_walkers());
  EXPECT_EQ(1, pars.get_accumulators());
  EXPECT_EQ(false, pars.shared_walk_and_accumulation_thread());
  EXPECT_EQ(false, pars.adjust_self_energy_for_double_counting());
  EXPECT_EQ(true, pars.store_configuration());
}

TEST(MciParametersTest, ReadAll) {
  dca::io::JSONReader reader;
  dca::phys::params::MciParameters pars;

  reader.open_file(DCA_SOURCE_DIR "/test/unit/phys/parameters/mci_parameters/input_read_all.json");
  pars.readWrite(reader);
  reader.close_file();

  EXPECT_EQ(42, pars.get_seed());
  EXPECT_EQ(40, pars.get_warm_up_sweeps());
  const auto expected_sweeps = std::vector<double>{1.5, 1.5, 4};
  EXPECT_EQ(expected_sweeps, pars.get_sweeps_per_measurement());
  const auto expected_meas = std::vector<int>{100, 150, 200};
  EXPECT_EQ(expected_meas, pars.get_measurements());
  EXPECT_EQ(dca::phys::ErrorComputationType::JACK_KNIFE, pars.get_error_computation_type());
  EXPECT_EQ(3, pars.get_walkers());
  EXPECT_EQ(5, pars.get_accumulators());
  EXPECT_EQ(true, pars.shared_walk_and_accumulation_thread());
  EXPECT_EQ(false, pars.store_configuration());
}

TEST(MciParametersTest, ReadPositiveIntegerSeed) {
  dca::io::JSONReader reader;
  dca::phys::params::MciParameters pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/mci_parameters/input_pos_int_seed.json");
  pars.readWrite(reader);
  reader.close_file();

  EXPECT_EQ(42, pars.get_seed());
}

TEST(MciParametersTest, ReadNegativeIntegerSeed) {
  dca::io::JSONReader reader;
  dca::phys::params::MciParameters pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/mci_parameters/input_neg_int_seed.json");
  pars.readWrite(reader);
  reader.close_file();

  EXPECT_EQ(-1, pars.get_seed());
}

TEST(MciParametersTest, ReadScalarForVector) {
  dca::io::JSONReader reader;
  dca::phys::params::MciParameters pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/mci_parameters/input_scalar_for_vec.json");
  pars.readWrite(reader);
  reader.close_file();

  EXPECT_EQ(std::vector<double>{4.}, pars.get_sweeps_per_measurement());
  EXPECT_EQ(std::vector<int>{200}, pars.get_measurements());
}

TEST(MciParametersTest, RandomSeed) {
  // The input file contains the seeding option "random" instead of a number.
  dca::io::JSONReader reader;
  dca::phys::params::MciParameters pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/mci_parameters/input_random_seed.json");

  std::vector<int> seeds;
  const int n_seeds = 5;

  for (int i = 0; i < n_seeds; i++) {
    pars.readWrite(reader);
    seeds.push_back(pars.get_seed());
  }

  reader.close_file();

  // Check that the seeds are in the expected interval.
  for (int i = 0; i < n_seeds; i++) {
    EXPECT_GE(seeds[i], 0);
    EXPECT_LE(seeds[i], std::numeric_limits<int>::max());
  }

  // Check that varying seeds are generated.
  std::sort(seeds.begin(), seeds.end());
  auto new_end = std::unique(seeds.begin(), seeds.end());
  EXPECT_GT(new_end - seeds.begin(), 1);
}

TEST(MciParametersTest, InvalidSeedingOption) {
  dca::io::JSONReader reader;
  dca::phys::params::MciParameters pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/mci_parameters/input_invalid_seeding_option.json");
  pars.readWrite(reader);
  reader.close_file();

  EXPECT_EQ(985456376, pars.get_seed());
}
