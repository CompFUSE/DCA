// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests function_parameters.hpp
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/phys/parameters/function_parameters.hpp"
#include "gtest/gtest.h"
#include "dca/io/json/json_reader.hpp"

TEST(FunctionParametersTest, DefaultValues) {
  int lattice_dim = 2;
  dca::phys::params::FunctionParameters pars(lattice_dim);

  EXPECT_EQ(std::vector<int>(lattice_dim, 8), pars.get_H_k_grid_size());
  EXPECT_EQ(128, pars.get_sp_time_intervals());
  EXPECT_EQ(256, pars.get_sp_fermionic_frequencies());
  EXPECT_EQ(32, pars.get_sp_bosonic_frequencies());
  EXPECT_EQ(std::vector<std::vector<int>>(0), pars.get_sp_cluster());
  EXPECT_EQ(-10, pars.get_min_real_frequency());
  EXPECT_EQ(10, pars.get_max_real_frequency());
  EXPECT_EQ(128, pars.get_number_of_real_frequencies());
  EXPECT_EQ(0.01, pars.get_real_frequencies_off_set());
  EXPECT_EQ(0, pars.get_tp_time_intervals());
  EXPECT_EQ(0, pars.get_tp_fermionic_frequencies());
  EXPECT_EQ(0, pars.get_tp_bosonic_frequencies());
  EXPECT_EQ(std::vector<std::vector<int>>(0), pars.get_tp_cluster());
}

TEST(FunctionParametersTest, ReadAll) {
  dca::io::JSONReader reader;
  int lattice_dim = 2;
  dca::phys::params::FunctionParameters pars(lattice_dim);

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/function_parameters/input_read_all.json");
  pars.readWrite(reader);
  reader.close_file();

  std::vector<std::vector<int>> sp_cluster_check{{4, 0}, {2, 4}};
  std::vector<std::vector<int>> tp_cluster_check{{4, 4}, {4, -4}};

  EXPECT_EQ(std::vector<int>(lattice_dim, 10), pars.get_H_k_grid_size());
  EXPECT_EQ(64, pars.get_sp_time_intervals());
  EXPECT_EQ(64, pars.get_sp_fermionic_frequencies());
  EXPECT_EQ(64, pars.get_sp_bosonic_frequencies());
  EXPECT_EQ(sp_cluster_check, pars.get_sp_cluster());
  EXPECT_EQ(-8, pars.get_min_real_frequency());
  EXPECT_EQ(8, pars.get_max_real_frequency());
  EXPECT_EQ(64, pars.get_number_of_real_frequencies());
  EXPECT_EQ(0.1, pars.get_real_frequencies_off_set());
  EXPECT_EQ(32, pars.get_tp_time_intervals());
  EXPECT_EQ(32, pars.get_tp_fermionic_frequencies());
  EXPECT_EQ(32, pars.get_tp_bosonic_frequencies());
  EXPECT_EQ(tp_cluster_check, pars.get_tp_cluster());
}
