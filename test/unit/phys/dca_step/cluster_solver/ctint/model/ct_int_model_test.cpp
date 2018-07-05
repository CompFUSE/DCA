// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
//

#include "gtest/gtest.h"

#include "dca/phys/dca_step/cluster_solver/ctint/model/general_interaction.hpp"

#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_rectangular.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/math/random/random.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/profiling/null_profiler.hpp"

const std::string input_dir =
    DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/ctint/inputs/";

TEST(CtIntModel, QuarticSquareLatticeTest) {
  using BaseLattice = dca::phys::models::square_lattice<dca::phys::domains::D4>;
  using LatticeType = dca::phys::models::GeneralInteraction<BaseLattice>;
  using Model = dca::phys::models::TightBindingModel<LatticeType>;
  using RngType = dca::math::random::StdRandomWrapper<std::ranlux48_base>;
  using Concurrency = dca::parallel::NoConcurrency;

  using ParametersType =
      dca::phys::params::Parameters<Concurrency, dca::parallel::NoThreading, dca::profiling::NullProfiler,
                                    Model, RngType, dca::phys::solver::CT_INT>;
  using Data = dca::phys::DcaData<ParametersType>;

  Concurrency concurrency(0, NULL);
  ParametersType parameters("", concurrency);

  parameters.read_input_and_broadcast<dca::io::JSONReader>(input_dir + "input.json");
  parameters.setHintFileName(input_dir + "H_int.txt");
  parameters.update_model();
  parameters.update_domains();

  Data data(parameters);
  data.initialize_H_0_and_H_i();
  EXPECT_EQ(2*4, data.H_interactions.integratedInteraction());
}
