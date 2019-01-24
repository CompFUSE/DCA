// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests output_parameters.hpp
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/phys/parameters/output_parameters.hpp"
#include "gtest/gtest.h"
#include "dca/io/json/json_reader.hpp"

TEST(OutputParametersTest, DefaultValues) {
  dca::phys::params::OutputParameters pars;

  EXPECT_EQ("./", pars.get_directory());
  EXPECT_EQ("HDF5", pars.get_output_format());
  EXPECT_EQ("", pars.get_directory_config_read());
  EXPECT_EQ("", pars.get_directory_config_write());
  EXPECT_EQ("dca.hdf5", pars.get_filename_dca());
  EXPECT_EQ("analysis.hdf5", pars.get_filename_analysis());
  EXPECT_EQ("ed.hdf5", pars.get_filename_ed());
  EXPECT_EQ("qmc.hdf5", pars.get_filename_qmc());
  EXPECT_EQ("profiling.json", pars.get_filename_profiling());
  EXPECT_FALSE(pars.dump_lattice_self_energy());
  EXPECT_FALSE(pars.dump_cluster_Greens_functions());
  EXPECT_FALSE(pars.dump_Gamma_lattice());
  EXPECT_FALSE(pars.dump_chi_0_lattice());
}

TEST(OutputParametersTest, ReadAll) {
  dca::io::JSONReader reader;
  dca::phys::params::OutputParameters pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/output_parameters/input_read_all.json");
  pars.readWrite(reader);
  reader.close_file();

  // HDF5 is the recommended output format. We use JSON in this test, since HDF5 is already the
  // default.
  EXPECT_EQ("./T=0.5", pars.get_directory());
  EXPECT_EQ("JSON", pars.get_output_format());
  EXPECT_EQ("configuration", pars.get_directory_config_read());
  EXPECT_EQ("configuration", pars.get_directory_config_write());
  EXPECT_EQ("dca.json", pars.get_filename_dca());
  EXPECT_EQ("analysis.json", pars.get_filename_analysis());
  EXPECT_EQ("ed.json", pars.get_filename_ed());
  EXPECT_EQ("qmc.json", pars.get_filename_qmc());
  EXPECT_EQ("profiling_run1.json", pars.get_filename_profiling());
  EXPECT_TRUE(pars.dump_lattice_self_energy());
  EXPECT_TRUE(pars.dump_cluster_Greens_functions());
  EXPECT_TRUE(pars.dump_Gamma_lattice());
  EXPECT_TRUE(pars.dump_chi_0_lattice());
}
