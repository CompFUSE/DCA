// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests analysis_parameters.hpp
//
// TODO: Add tests for get_buffer_size, pack, unpack, and writing.

#include "dca/phys/parameters/analysis_parameters.hpp"
#include "gtest/gtest.h"
#include "dca/io/json/json_reader.hpp"

TEST(AnalysisParametersTest, DefaultValues) {
  dca::phys::params::AnalysisParameters pars;

  EXPECT_TRUE(pars.symmetrize_Gamma());
  EXPECT_EQ(0.5, pars.get_Gamma_deconvolution_cut_off());
  EXPECT_FALSE(pars.project_onto_crystal_harmonics());
  EXPECT_EQ(1.5, pars.get_projection_cut_off_radius());
}

TEST(AnalysisParametersTest, ReadAll) {
  dca::io::JSONReader reader;
  dca::phys::params::AnalysisParameters pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/analysis_parameters/input_read_all.json");
  pars.readWrite(reader);
  reader.close_file();

  std::vector<double> momentum_transfer_input_check{3.14, -1.57};

  EXPECT_FALSE(pars.symmetrize_Gamma());
  EXPECT_EQ(0.1, pars.get_Gamma_deconvolution_cut_off());
  EXPECT_TRUE(pars.project_onto_crystal_harmonics());
  EXPECT_EQ(3.0, pars.get_projection_cut_off_radius());
}
