// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests domains_parameters.hpp
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/phys/parameters/domains_parameters.hpp"
#include <vector>
#include "gtest/gtest.h"
#include "dca/io/json/json_reader.hpp"

TEST(DomainsParametersTest, DefaultValues) {
  const int dimension = 2;
  dca::phys::params::DomainsParameters pars(dimension);

  const std::vector<std::vector<int>> cluster_check{{1, 0}, {0, 1}};

  EXPECT_EQ(cluster_check, pars.get_cluster());
  EXPECT_EQ(1, pars.get_sp_superlattice_size());
  EXPECT_EQ(cluster_check, pars.get_sp_host());
  EXPECT_EQ(1, pars.get_tp_superlattice_size());
  EXPECT_EQ(cluster_check, pars.get_tp_host());
  EXPECT_EQ(128, pars.get_sp_time_intervals());
  EXPECT_EQ(1, pars.get_time_intervals_for_time_measurements());
  EXPECT_EQ(256, pars.get_sp_fermionic_frequencies());
  EXPECT_EQ(0, pars.get_hts_bosonic_frequencies());
  EXPECT_EQ(1, pars.get_four_point_fermionic_frequencies());
  EXPECT_EQ(-10., pars.get_min_real_frequency());
  EXPECT_EQ(10., pars.get_max_real_frequency());
  EXPECT_EQ(3, pars.get_real_frequencies());
  EXPECT_EQ(0.01, pars.get_imaginary_damping());
}

TEST(DomainsParametersTest, ReadAll) {
  dca::io::JSONReader reader;
  const int dimension = 2;
  dca::phys::params::DomainsParameters pars(dimension);

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/domains_parameters/input_read_all.json");
  pars.readWrite(reader);
  reader.close_file();

  const std::vector<std::vector<int>> cluster_check{{2, 0}, {0, 2}};
  const std::vector<std::vector<int>> sp_host_check{{20, 0}, {0, 20}};
  const std::vector<std::vector<int>> tp_host_check{{12, 0}, {0, 12}};

  EXPECT_EQ(cluster_check, pars.get_cluster());
  EXPECT_EQ(10, pars.get_sp_superlattice_size());
  EXPECT_EQ(sp_host_check, pars.get_sp_host());
  EXPECT_EQ(6, pars.get_tp_superlattice_size());
  EXPECT_EQ(tp_host_check, pars.get_tp_host());
  EXPECT_EQ(64, pars.get_sp_time_intervals());
  EXPECT_EQ(128, pars.get_time_intervals_for_time_measurements());
  EXPECT_EQ(512, pars.get_sp_fermionic_frequencies());
  EXPECT_EQ(32, pars.get_hts_bosonic_frequencies());
  EXPECT_EQ(16, pars.get_four_point_fermionic_frequencies());
  EXPECT_EQ(-8., pars.get_min_real_frequency());
  EXPECT_EQ(8., pars.get_max_real_frequency());
  EXPECT_EQ(128, pars.get_real_frequencies());
  EXPECT_EQ(0.001, pars.get_imaginary_damping());
}

TEST(DomainsParametersTest, ReadInvalid) {
  dca::io::JSONReader reader;
  const int dimension = 2;
  dca::phys::params::DomainsParameters pars(dimension);

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/domains_parameters/input_read_invalid.json");
  pars.readWrite(reader);
  reader.close_file();

  EXPECT_EQ(4, pars.get_sp_superlattice_size());
  EXPECT_EQ(2, pars.get_tp_superlattice_size());
}

TEST(DomainsParametersTest, Setters) {
  const int dimension = 2;
  dca::phys::params::DomainsParameters pars(dimension);

  // cluster
  const std::vector<std::vector<int>> cluster{{2, 0}, {0, 2}};
  pars.set_cluster(cluster);
  EXPECT_EQ(cluster, pars.get_cluster());

  const std::vector<std::vector<int>> cluster_wrong_dim{{2, 0, 0}, {0, 2, 0}, {0, 0, 2}};
  EXPECT_THROW(pars.set_cluster(cluster_wrong_dim), std::invalid_argument);

  // sp_superlattice_size
  const int sp_superlattice_size = 10;
  pars.set_sp_superlattice_size(sp_superlattice_size);
  EXPECT_EQ(sp_superlattice_size, pars.get_sp_superlattice_size());

  const int sp_superlattice_size_odd = 3;
  pars.set_sp_superlattice_size(sp_superlattice_size_odd);
  EXPECT_EQ(sp_superlattice_size_odd - 1, pars.get_sp_superlattice_size());

  // tp_superlattice_size
  const int tp_superlattice_size = 10;
  pars.set_tp_superlattice_size(tp_superlattice_size);
  EXPECT_EQ(tp_superlattice_size, pars.get_tp_superlattice_size());

  const int tp_superlattice_size_odd = 3;
  pars.set_tp_superlattice_size(tp_superlattice_size_odd);
  EXPECT_EQ(tp_superlattice_size_odd - 1, pars.get_tp_superlattice_size());
}
