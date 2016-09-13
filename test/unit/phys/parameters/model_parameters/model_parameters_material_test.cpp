// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests the material models specialization of the ModelParameters class.
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/phys/parameters/model_parameters.hpp"
#include "gtest/gtest.h"
#include "comp_library/IO_library/JSON/JSON.hpp"
#include "phys_library/domains/cluster/symmetries/point_groups/2D/2D_square.h"

using PointGroup = D4;

TEST(ModelParametersMaterialTest, DefaultValues) {
  dca::phys::params::ModelParameters<
      dca::phys::models::TightBindingModel<material_lattice<NiO_symmetric, PointGroup>>>
      pars;

  EXPECT_EQ("t_ij_file_name", pars.get_t_ij_file_name());
  EXPECT_EQ("U_ij_file_name", pars.get_U_ij_file_name());
}

TEST(ModelParametersMaterialTest, ReadAll) {
  IO::reader<IO::JSON> reader;
  dca::phys::params::ModelParameters<
      dca::phys::models::TightBindingModel<material_lattice<NiO_symmetric, PointGroup>>>
      pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/model_parameters/input_read_all_material.json");
  pars.readWrite(reader);
  reader.close_file();

  EXPECT_EQ("t_ij.txt", pars.get_t_ij_file_name());
  EXPECT_EQ("U_ij.txt", pars.get_U_ij_file_name());
}

TEST(ModelParametersMaterialTest, Setter) {
  dca::phys::params::ModelParameters<
      dca::phys::models::TightBindingModel<material_lattice<NiO_symmetric, PointGroup>>>
      pars;

  pars.set_t_ij_file_name("t_ij_NiO.txt");
  pars.set_U_ij_file_name("U_ij_NiO.txt");

  EXPECT_EQ("t_ij_NiO.txt", pars.get_t_ij_file_name());
  EXPECT_EQ("U_ij_NiO.txt", pars.get_U_ij_file_name());
}
