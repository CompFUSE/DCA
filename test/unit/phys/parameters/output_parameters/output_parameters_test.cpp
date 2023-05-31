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


#include "dca/config/haves_defines.hpp"
#include "dca/phys/parameters/output_parameters.hpp"
#include "dca/testing/gtest_h_w_warning_blocking.h"
#include "dca/io/json/json_reader.hpp"

TEST(OutputParametersTest, DefaultValues) {
  dca::phys::params::OutputParameters pars;

  EXPECT_EQ("./", pars.get_directory());
#ifdef DCA_HAVE_ADIOS2
  EXPECT_EQ("ADIOS2", pars.get_output_format());
  EXPECT_EQ("ADIOS2", pars.get_g4_output_format());
#else
  EXPECT_EQ("HDF5", pars.get_output_format());
  EXPECT_EQ("", pars.get_g4_output_format());
#endif
  EXPECT_FALSE(pars.autoresume());
  EXPECT_EQ("", pars.get_directory_config_read());
  EXPECT_EQ("", pars.get_directory_config_write());
  EXPECT_EQ("dca.bp", pars.get_filename_dca());
  EXPECT_EQ("sofqomega.bp", pars.get_filename_analysis());
  EXPECT_EQ("ed.hdf5", pars.get_filename_ed());
  EXPECT_EQ("qmc.hdf5", pars.get_filename_qmc());
  EXPECT_EQ("profiling.json", pars.get_filename_profiling());
  EXPECT_FALSE(pars.dump_lattice_self_energy());
  EXPECT_FALSE(pars.dump_cluster_Greens_functions());
  EXPECT_FALSE(pars.dump_Gamma_lattice());
  EXPECT_FALSE(pars.dump_chi_0_lattice());
  EXPECT_FALSE(pars.dump_every_iteration());

  // // And nothing except required input should get updated by reading the minimal input
  // dca::io::JSONReader reader;
  // reader.open_file(DCA_SOURCE_DIR
  //               "/test/unit/phys/parameters/output_parameters/input_check_defaults.json");
  // pars.readWrite(reader);
  // reader.close_file();
  // EXPECT_EQ("./", pars.get_directory());
  // EXPECT_EQ("ADIOS2", pars.get_output_format());
  // EXPECT_EQ("ADIOS2", pars.get_g4_output_format()),
  // EXPECT_EQ(false, pars.autoresume());
  // EXPECT_EQ("", pars.get_directory_config_read());
  // EXPECT_EQ("", pars.get_directory_config_write());
  // EXPECT_EQ("dca.bp", pars.get_filename_dca());
  // EXPECT_EQ("sofqomega.bp", pars.get_filename_analysis());
  // EXPECT_EQ("ed.hdf5", pars.get_filename_ed());
  // EXPECT_EQ("qmc.hdf5", pars.get_filename_qmc());
  // EXPECT_EQ("profiling.json", pars.get_filename_profiling());
  // EXPECT_FALSE(pars.dump_lattice_self_energy());
  // EXPECT_FALSE(pars.dump_cluster_Greens_functions());
  // EXPECT_FALSE(pars.dump_Gamma_lattice());
  // EXPECT_FALSE(pars.dump_chi_0_lattice());
  // EXPECT_FALSE(pars.dump_every_iteration());
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
  EXPECT_EQ("HDF5", pars.get_output_format());
  EXPECT_EQ(true, pars.autoresume());
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
