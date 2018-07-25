// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests dca_parameters.hpp
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/phys/parameters/dca_parameters.hpp"
#include "gtest/gtest.h"
#include "dca/io/json/json_reader.hpp"

class DcaParametersTest : public ::testing::Test {
protected:
  dca::phys::params::DcaParameters pars_;
  dca::io::JSONReader reader_;
};

TEST_F(DcaParametersTest, DefaultValues) {
  EXPECT_EQ("zero", pars_.get_initial_self_energy());
  EXPECT_EQ(1, pars_.get_dca_iterations());
  EXPECT_EQ(0., pars_.get_dca_accuracy());
  EXPECT_EQ(1., pars_.get_self_energy_mixing_factor());
  EXPECT_EQ(std::vector<int>{0}, pars_.get_interacting_orbitals());
  EXPECT_FALSE(pars_.do_finite_size_qmc());
  EXPECT_EQ(true, pars_.do_simple_q_points_summation());
  EXPECT_EQ(0, pars_.get_k_mesh_recursion());
  EXPECT_EQ(0, pars_.get_coarsegraining_periods());
  EXPECT_EQ(1, pars_.get_quadrature_rule());
  EXPECT_EQ(1, pars_.get_coarsegraining_threads());
  EXPECT_EQ(0, pars_.get_tail_frequencies());
  EXPECT_FALSE(pars_.do_dca_plus());
  EXPECT_EQ(16, pars_.get_deconvolution_iterations());
  EXPECT_EQ(1.e-3, pars_.get_deconvolution_tolerance());
  EXPECT_FALSE(pars_.hts_approximation());
  EXPECT_EQ(1, pars_.get_hts_threads());
}

TEST_F(DcaParametersTest, ReadAll) {
  reader_.open_file(DCA_SOURCE_DIR "/test/unit/phys/parameters/dca_parameters/input_read_all.json");
  pars_.readWrite(reader_);
  reader_.close_file();

  std::vector<int> interacting_orbitals_check{0, 1, 2};

  EXPECT_EQ("./T=0.5/dca.hdf5", pars_.get_initial_self_energy());
  EXPECT_EQ(3, pars_.get_dca_iterations());
  EXPECT_EQ(1.e-3, pars_.get_dca_accuracy());
  EXPECT_EQ(0.5, pars_.get_self_energy_mixing_factor());
  EXPECT_EQ(interacting_orbitals_check, pars_.get_interacting_orbitals());
  EXPECT_FALSE(pars_.do_finite_size_qmc());
  EXPECT_FALSE(pars_.do_simple_q_points_summation());
  EXPECT_EQ(3, pars_.get_k_mesh_recursion());
  EXPECT_EQ(2, pars_.get_coarsegraining_periods());
  EXPECT_EQ(2, pars_.get_quadrature_rule());
  EXPECT_EQ(8, pars_.get_coarsegraining_threads());
  EXPECT_EQ(10, pars_.get_tail_frequencies());
  EXPECT_TRUE(pars_.do_dca_plus());
  EXPECT_EQ(32, pars_.get_deconvolution_iterations());
  EXPECT_EQ(1.e-4, pars_.get_deconvolution_tolerance());
  EXPECT_TRUE(pars_.hts_approximation());
  EXPECT_EQ(8, pars_.get_hts_threads());
}

// Separate test for reading the "do-finite-size-QMC" parameter, since we cannot set it to true at
// the same time as "do-DCA+" is set to true.
TEST_F(DcaParametersTest, ReadFiniteSizeQMC) {
  reader_.open_file(DCA_SOURCE_DIR
                    "/test/unit/phys/parameters/dca_parameters/finite_size_qmc.json");
  pars_.readWrite(reader_);
  reader_.close_file();

  EXPECT_TRUE(pars_.do_finite_size_qmc());
}

TEST_F(DcaParametersTest, ConsistencyCheck) {
  reader_.open_file(DCA_SOURCE_DIR "/test/unit/phys/parameters/dca_parameters/not_consistent.json");
  EXPECT_THROW(pars_.readWrite(reader_), std::logic_error);
  reader_.close_file();
}
