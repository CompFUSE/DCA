// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests ed_solver_parameters.hpp
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/phys/parameters/ed_solver_parameters.hpp"
#include "gtest/gtest.h"
#include "comp_library/IO_library/JSON/JSON.hpp"

TEST(EdSolverParametersTest, DefaultValues) {
  dca::phys::params::EdSolverParameters pars;

  EXPECT_EQ(1.e-6, pars.get_eigenvalue_cut_off());
  EXPECT_EQ("default", pars.get_ed_method());
  EXPECT_EQ(0, pars.get_occupation());
  EXPECT_EQ(0, pars.get_magnetization());
  EXPECT_FALSE(pars.check_orthogonality_of_states());
}

TEST(EdSolverParametersTest, ReadAll) {
  IO::reader<IO::JSON> reader;
  dca::phys::params::EdSolverParameters pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/ed_solver_parameters/input_read_all.json");
  pars.readWrite(reader);
  reader.close_file();

  EXPECT_EQ(1.e-8, pars.get_eigenvalue_cut_off());
  EXPECT_EQ("block-diagonal", pars.get_ed_method());
  EXPECT_EQ(1, pars.get_occupation());
  EXPECT_EQ(-1, pars.get_magnetization());
  EXPECT_TRUE(pars.check_orthogonality_of_states());
}
