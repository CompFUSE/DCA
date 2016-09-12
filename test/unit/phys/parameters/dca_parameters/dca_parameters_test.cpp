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
#include "comp_library/IO_library/JSON/JSON.hpp"

TEST(DcaParametersTest, DefaultValues) {
  int lattice_dim = 2;
  dca::phys::params::DcaParameters pars(lattice_dim);

  std::vector<std::vector<int>> DCA_cluster_check{{0, 0}, {0, 0}};

  EXPECT_EQ(false, pars.do_DCA_plus());
  EXPECT_EQ(std::vector<int>(0), pars.get_interacting_bands());
  EXPECT_EQ(1, pars.get_DCA_iterations());
  EXPECT_EQ(1.e-5, pars.get_DCA_accuracy());
  EXPECT_EQ(1., pars.get_DCA_mixing_factor());
  EXPECT_EQ(DCA_cluster_check, pars.get_DCA_cluster());
  EXPECT_EQ(3, pars.get_k_mesh_refinement());
  EXPECT_EQ(0, pars.get_number_of_periods());
  EXPECT_EQ(3, pars.get_quadrature_rule());
  EXPECT_EQ(true, pars.precompute_Hamiltonian());
  EXPECT_EQ(1, pars.get_nr_coarsegraining_threads());
  EXPECT_EQ(0, pars.get_number_of_tail_frequencies());
  EXPECT_EQ(1.e-3, pars.get_phi_k_integration_accuracy());
  EXPECT_EQ(false, pars.print_phi_k());
  EXPECT_EQ("wannier-interpolation", pars.get_interpolation_method());
  EXPECT_EQ(false, pars.use_HTS_approximation());
  EXPECT_EQ(1.e-2, pars.get_deconvolution_tolerance());
  EXPECT_EQ(16, pars.get_max_deconvolution_iterations());
}

TEST(DcaParametersTest, ReadAll) {
  IO::reader<IO::JSON> reader;
  int lattice_dim = 2;
  dca::phys::params::DcaParameters pars(lattice_dim);

  reader.open_file(DCA_SOURCE_DIR "/test/unit/phys/parameters/dca_parameters/input_read_all.json");
  pars.readWrite(reader);
  reader.close_file();

  std::vector<int> interacting_bands_check{0, 1};
  std::vector<std::vector<int>> DCA_cluster_check{{4, 0}, {0, 4}};

  EXPECT_EQ(true, pars.do_DCA_plus());
  EXPECT_EQ(interacting_bands_check, pars.get_interacting_bands());
  EXPECT_EQ(3, pars.get_DCA_iterations());
  EXPECT_EQ(1.e-6, pars.get_DCA_accuracy());
  EXPECT_EQ(0.5, pars.get_DCA_mixing_factor());
  EXPECT_EQ(DCA_cluster_check, pars.get_DCA_cluster());
  EXPECT_EQ(4, pars.get_k_mesh_refinement());
  EXPECT_EQ(2, pars.get_number_of_periods());
  EXPECT_EQ(2, pars.get_quadrature_rule());
  EXPECT_EQ(false, pars.precompute_Hamiltonian());
  EXPECT_EQ(8, pars.get_nr_coarsegraining_threads());
  EXPECT_EQ(10, pars.get_number_of_tail_frequencies());
  EXPECT_EQ(1.e-5, pars.get_phi_k_integration_accuracy());
  EXPECT_EQ(true, pars.print_phi_k());
  EXPECT_EQ("some-method", pars.get_interpolation_method());
  EXPECT_EQ(true, pars.use_HTS_approximation());
  EXPECT_EQ(1.e-3, pars.get_deconvolution_tolerance());
  EXPECT_EQ(32, pars.get_max_deconvolution_iterations());

  EXPECT_TRUE(pars.is_an_interacting_band(1));
  EXPECT_FALSE(pars.is_an_interacting_band(3));
}
