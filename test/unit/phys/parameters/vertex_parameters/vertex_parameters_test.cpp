// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests vertex_parameters.hpp
//
// TODO: Add tests for get_buffer_size, pack, unpack, writing, and the computation of the 'exact'
//       q-channel vector and its index.

#include "dca/phys/parameters/vertex_parameters.hpp"
#include "gtest/gtest.h"
#include "dca/io/json/json_reader.hpp"

TEST(VertexParametersTest, DefaultValues) {
  dca::phys::params::VertexParameters<2> pars;

  EXPECT_EQ(NONE, pars.get_vertex_measurement_type());
  EXPECT_EQ(std::vector<double>(2, 0.), pars.get_q_channel_vec_input());
  EXPECT_EQ(0, pars.get_w_channel());
  EXPECT_EQ(0.5, pars.get_singular_value_cut_off());
  EXPECT_EQ(128, pars.get_singular_value_index_cut_off());
  EXPECT_FALSE(pars.do_diagonalization_on_folded_Gamma_chi_0());
  EXPECT_EQ(10., pars.get_BSE_cut_off_radius());
  EXPECT_TRUE(pars.do_deconvolution_of_Gamma());
  EXPECT_TRUE(pars.do_symmetrization_of_Gamma());
  EXPECT_FALSE(pars.compute_chi_0());
  EXPECT_FALSE(pars.compute_chi());
  EXPECT_TRUE(pars.compute_eigenvalues());
  EXPECT_FALSE(pars.compute_P_q_cluster());
  EXPECT_FALSE(pars.compute_P_q_lattice());
}

TEST(VertexParametersTest, ReadAll) {
  dca::io::JSONReader reader;
  dca::phys::params::VertexParameters<2> pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/vertex_parameters/input_read_all.json");
  pars.readWrite(reader);
  reader.close_file();

  std::vector<double> q_channel_vec_check{3.14, -1.57};

  EXPECT_EQ(PARTICLE_PARTICLE_SUPERCONDUCTING, pars.get_vertex_measurement_type());
  EXPECT_EQ(q_channel_vec_check, pars.get_q_channel_vec_input());
  EXPECT_EQ(2, pars.get_w_channel());
  EXPECT_EQ(1.3, pars.get_singular_value_cut_off());
  EXPECT_EQ(99, pars.get_singular_value_index_cut_off());
  EXPECT_TRUE(pars.do_diagonalization_on_folded_Gamma_chi_0());
  EXPECT_EQ(7.5, pars.get_BSE_cut_off_radius());
  EXPECT_FALSE(pars.do_deconvolution_of_Gamma());
  EXPECT_FALSE(pars.do_symmetrization_of_Gamma());
  EXPECT_TRUE(pars.compute_chi_0());
  EXPECT_TRUE(pars.compute_chi());
  EXPECT_FALSE(pars.compute_eigenvalues());
  EXPECT_TRUE(pars.compute_P_q_cluster());
  EXPECT_TRUE(pars.compute_P_q_lattice());
}
