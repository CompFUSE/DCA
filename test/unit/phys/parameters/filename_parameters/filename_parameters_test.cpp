// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests filename_parameters.hpp
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/phys/parameters/filename_parameters.hpp"
#include "gtest/gtest.h"
#include "dca/io/json/json_reader.hpp"

TEST(FilenameParametersTest, DefaultValues) {
  dca::phys::params::FilenameParameters pars;

  EXPECT_EQ("./", pars.get_directory());
  EXPECT_EQ("JSON", pars.get_output_format());
  EXPECT_EQ(".json", pars.get_file_extension());
  EXPECT_EQ("output.json", pars.get_output_file_name());
  EXPECT_EQ("prof_data.txt", pars.get_profiling_file_name());
  EXPECT_EQ("data_spectrum.json", pars.get_spectrum_file_name());
  EXPECT_EQ("data_susceptibilities.json", pars.get_susceptibilities_file_name());
  EXPECT_EQ("vertex_filename.py", pars.get_vertex_file_name());
  EXPECT_EQ("output_ED.json", pars.get_ED_output_file_name());
  EXPECT_EQ("output_CPE.json", pars.get_CPE_output_file_name());
  EXPECT_EQ("output_QMC.json", pars.get_QMC_output_file_name());
  EXPECT_FALSE(pars.dump_lattice_Self_energy());
  EXPECT_FALSE(pars.dump_cluster_Greens_functions());
}

TEST(FilenameParametersTest, ReadAll) {
  dca::io::JSONReader reader;
  dca::phys::params::FilenameParameters pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/filename_parameters/input_read_all.json");
  pars.readWrite(reader);
  reader.close_file();

  EXPECT_EQ("./T=0.5", pars.get_directory());
  EXPECT_EQ("HDF5", pars.get_output_format());
  EXPECT_EQ(".hdf5", pars.get_file_extension());
  EXPECT_EQ("output_run1.hdf5", pars.get_output_file_name());
  EXPECT_EQ("profiling_run1.txt", pars.get_profiling_file_name());
  EXPECT_EQ("spectrum_run1.json", pars.get_spectrum_file_name());
  EXPECT_EQ("susceptibilities_run1.json", pars.get_susceptibilities_file_name());
  EXPECT_EQ("plot_vertex.py", pars.get_vertex_file_name());
  EXPECT_EQ("output_ed_run1.json", pars.get_ED_output_file_name());
  EXPECT_EQ("output_cpe_run1.json", pars.get_CPE_output_file_name());
  EXPECT_EQ("output_qmc_run1.json", pars.get_QMC_output_file_name());
  EXPECT_TRUE(pars.dump_lattice_Self_energy());
  EXPECT_TRUE(pars.dump_cluster_Greens_functions());
}
