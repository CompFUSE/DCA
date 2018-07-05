// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file tests the SS-CT-HYB specialization of the McSolverParameters class.
//
// TODO: Add tests for get_buffer_size, pack, unpack and writing.

#include "dca/phys/parameters/mc_solver_parameters.hpp"
#include "gtest/gtest.h"
#include "dca/io/json/json_reader.hpp"

TEST(McSolverParametersSsCtHybTest, DefaultValues) {
  dca::phys::params::McSolverParameters<dca::phys::solver::SS_CT_HYB> pars;

  EXPECT_EQ(0, pars.get_self_energy_tail_cutoff());
  EXPECT_EQ(0.5, pars.get_steps_per_sweep());
  EXPECT_EQ(0.5, pars.get_shifts_per_sweep());
}

TEST(McSolverParametersSsCtHybTest, ReadAll) {
  dca::io::JSONReader reader;
  dca::phys::params::McSolverParameters<dca::phys::solver::SS_CT_HYB> pars;

  reader.open_file(DCA_SOURCE_DIR
                   "/test/unit/phys/parameters/mc_solver_parameters/input_read_all_ss_ct_hyb.json");
  pars.readWrite(reader);
  reader.close_file();

  EXPECT_EQ(10, pars.get_self_energy_tail_cutoff());
  EXPECT_EQ(0.6, pars.get_steps_per_sweep());
  EXPECT_EQ(0.4, pars.get_shifts_per_sweep());
}
