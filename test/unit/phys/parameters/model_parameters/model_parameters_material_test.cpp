// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests the material models specialization of the ModelParameters class.
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/phys/parameters/model_parameters.hpp"
#include "gtest/gtest.h"
#include "dca/io/json/json_reader.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"

using PointGroup = dca::phys::domains::D4;

TEST(ModelParametersMaterialTest, DefaultValues) {
  dca::phys::params::ModelParameters<dca::phys::models::TightBindingModel<
      dca::phys::models::material_lattice<dca::phys::models::NiO_symmetric, PointGroup>>>
      pars;

  EXPECT_EQ("t_ij.txt", pars.get_t_ij_file_name());
  EXPECT_EQ("U_ij.txt", pars.get_U_ij_file_name());
}

TEST(ModelParametersMaterialTest, ReadAll) {
  dca::io::JSONReader reader;
  dca::phys::params::ModelParameters<dca::phys::models::TightBindingModel<
      dca::phys::models::material_lattice<dca::phys::models::NiO_symmetric, PointGroup>>>
      pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/model_parameters/input_read_all_material.json");
  pars.readWrite(reader);
  reader.close_file();

  EXPECT_EQ("NiO_t_ij.txt", pars.get_t_ij_file_name());
  EXPECT_EQ("NiO_U_ij.txt", pars.get_U_ij_file_name());
}

TEST(ModelParametersMaterialTest, Setter) {
  dca::phys::params::ModelParameters<dca::phys::models::TightBindingModel<
      dca::phys::models::material_lattice<dca::phys::models::NiO_symmetric, PointGroup>>>
      pars;

  pars.set_t_ij_file_name("NiO_t_ij.txt");
  pars.set_U_ij_file_name("NiO_U_ij.txt");

  EXPECT_EQ("NiO_t_ij.txt", pars.get_t_ij_file_name());
  EXPECT_EQ("NiO_U_ij.txt", pars.get_U_ij_file_name());
}
