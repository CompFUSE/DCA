// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests ed_solver_parameters.hpp
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/phys/parameters/ed_solver_parameters.hpp"
#include "gtest/gtest.h"
#include "dca/io/json/json_reader.hpp"

TEST(EdSolverParametersTest, DefaultValues) {
  dca::phys::params::EdSolverParameters pars;
  EXPECT_EQ(1.e-6, pars.get_eigenvalue_cut_off());
}

TEST(EdSolverParametersTest, ReadAll) {
  dca::io::JSONReader reader;
  dca::phys::params::EdSolverParameters pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/ed_solver_parameters/input_read_all.json");
  pars.readWrite(reader);
  reader.close_file();

  EXPECT_EQ(1.e-4, pars.get_eigenvalue_cut_off());
}
