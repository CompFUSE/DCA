// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests dca_parameters.hpp
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/phys/parameters/dca_parameters.hpp"
#include "gtest/gtest.h"
#include "dca/io/json/json_reader.hpp"

TEST(DcaParametersTest, DefaultValues) {
  dca::phys::params::DcaParameters pars;

  EXPECT_EQ("zero", pars.get_initial_self_energy());
  EXPECT_EQ(1, pars.get_dca_iterations());
  EXPECT_EQ(0., pars.get_dca_accuracy());
  EXPECT_EQ(1., pars.get_self_energy_mixing_factor());
  EXPECT_EQ(std::vector<int>{0}, pars.get_interacting_orbitals());
  EXPECT_EQ(std::vector<std::vector<int>>{}, pars.get_cluster());
  EXPECT_EQ(std::vector<std::vector<int>>{}, pars.get_sp_host());
  EXPECT_EQ(std::vector<std::vector<int>>{}, pars.get_tp_host());
  EXPECT_EQ(128, pars.get_sp_time_intervals());
  EXPECT_EQ(256, pars.get_sp_fermionic_frequencies());
  EXPECT_EQ(0, pars.get_k_mesh_recursion());
  EXPECT_EQ(0, pars.get_coarsegraining_periods());
  EXPECT_EQ(1, pars.get_quadrature_rule());
  EXPECT_EQ(1, pars.get_coarsegraining_threads());
  EXPECT_EQ(0, pars.get_tail_frequencies());
  EXPECT_FALSE(pars.hts_approximation());
  EXPECT_EQ(32, pars.get_hts_bosonic_frequencies());
  EXPECT_EQ(1, pars.get_hts_threads());
  EXPECT_FALSE(pars.do_dca_plus());
  EXPECT_EQ(16, pars.get_deconvolution_iterations());
  EXPECT_EQ(1.e-3, pars.get_deconvolution_tolerance());
}

TEST(DcaParametersTest, ReadAll) {
  dca::io::JSONReader reader;
  dca::phys::params::DcaParameters pars;

  reader.open_file(DCA_SOURCE_DIR "/test/unit/phys/parameters/dca_parameters/input_read_all.json");
  pars.readWrite(reader);
  reader.close_file();

  std::vector<int> interacting_orbitals_check{0, 1, 2};
  std::vector<std::vector<int>> cluster_check{{2, 0}, {0, 2}};
  std::vector<std::vector<int>> sp_host_check{{20, 20}, {20, -20}};
  std::vector<std::vector<int>> tp_host_check{{8, 8}, {8, -8}};

  EXPECT_EQ("./T=0.5/dca.hdf5", pars.get_initial_self_energy());
  EXPECT_EQ(3, pars.get_dca_iterations());
  EXPECT_EQ(1.e-3, pars.get_dca_accuracy());
  EXPECT_EQ(0.5, pars.get_self_energy_mixing_factor());
  EXPECT_EQ(interacting_orbitals_check, pars.get_interacting_orbitals());
  EXPECT_EQ(cluster_check, pars.get_cluster());
  EXPECT_EQ(sp_host_check, pars.get_sp_host());
  EXPECT_EQ(tp_host_check, pars.get_tp_host());
  EXPECT_EQ(64, pars.get_sp_time_intervals());
  EXPECT_EQ(128, pars.get_sp_fermionic_frequencies());
  EXPECT_EQ(3, pars.get_k_mesh_recursion());
  EXPECT_EQ(2, pars.get_coarsegraining_periods());
  EXPECT_EQ(2, pars.get_quadrature_rule());
  EXPECT_EQ(8, pars.get_coarsegraining_threads());
  EXPECT_EQ(10, pars.get_tail_frequencies());
  EXPECT_TRUE(pars.hts_approximation());
  EXPECT_EQ(16, pars.get_hts_bosonic_frequencies());
  EXPECT_EQ(8, pars.get_hts_threads());
  EXPECT_TRUE(pars.do_dca_plus());
  EXPECT_EQ(32, pars.get_deconvolution_iterations());
  EXPECT_EQ(1.e-4, pars.get_deconvolution_tolerance());
}
