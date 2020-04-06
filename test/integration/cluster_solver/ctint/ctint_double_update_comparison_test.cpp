// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// No-change test for CT-INT.
// Bilayer lattice with two band and two sites.

#include <iostream>
#include <string>

#include "gtest/gtest.h"

#include "dca/phys/dca_step/cluster_solver/ctint/details/solver_methods.hpp"
#include "dca/function/function.hpp"
#include "dca/function/util/difference.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "test/unit/phys/dca_step/cluster_solver/stub_rng.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/parallel//no_concurrency/no_concurrency.hpp"
#include "test/unit/phys/dca_step/cluster_solver/ctint/walker/walker_wrapper_submatrix.hpp"
#include "test/unit/phys/dca_step/cluster_solver/ctint/walker/walker_wrapper.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/fe_as_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"

const std::string input_dir = DCA_SOURCE_DIR "/test/integration/cluster_solver/ctint/";

template <class WalkerType, class G0, class Parameters, class Data>
void initializeWalkerStatic(const G0& g0, const Parameters& parameters, const Data& data) {
  WalkerType::setDMatrixBuilder(g0);
  WalkerType::setDMatrixAlpha(parameters.getAlphas(), parameters.adjustAlphaDd());
  WalkerType::setInteractionVertices(data, parameters);
}

TEST(CtintDoubleUpdateComparisonTest, Self_Energy) {
  using RngType = dca::testing::StubRng;
  using RealRng = dca::math::random::StdRandomWrapper<std::mt19937_64>;
  using Lattice = dca::phys::models::FeAsLattice<dca::phys::domains::D4>;
  using Model = dca::phys::models::TightBindingModel<Lattice>;
  using Threading = dca::parallel::NoThreading;
  using Concurrency = dca::parallel::NoConcurrency;
  using Parameters =
      dca::phys::params::Parameters<Concurrency, Threading, dca::profiling::NullProfiler, Model,
                                    RngType, dca::phys::solver::CT_INT>;
  using Data = dca::phys::DcaData<Parameters>;

  using Walker = testing::phys::solver::ctint::WalkerWrapper<Parameters, double>;
  using WalkerSubmatrix =
      testing::phys::solver::ctint::WalkerWrapperSubmatrix<Parameters, dca::linalg::CPU, double>;

  Concurrency concurrency(0, nullptr);
  dca::util::GitVersion::print();
  dca::util::Modules::print();

  Parameters parameters(dca::util::GitVersion::string(), concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(
      input_dir + "/double_insertion_comparison_input.json");
  parameters.update_model();
  parameters.update_domains();

  Data data(parameters);
  data.initialize();

  dca::phys::solver::ctint::G0Interpolation<dca::linalg::CPU, double> g0(
      dca::phys::solver::ctint::details::shrinkG0(data.G0_r_t));

  initializeWalkerStatic<Walker>(g0, parameters, data);
  initializeWalkerStatic<WalkerSubmatrix>(g0, parameters, data);

  RealRng rng(0, 1);
  std::vector<double> rng_vals(10000);
  for (auto& x : rng_vals)
    x = rng();
  RngType rng1(rng_vals), rng2(rng_vals);

  Walker walker1(parameters, rng1);

  parameters.setMaxSubmatrixSize(16);
  WalkerSubmatrix walker2(parameters, rng2);

  EXPECT_NEAR(walker1.get_MC_log_weight(), walker2.get_MC_log_weight(), 5e-7);

  for (int i = 0; i < 128; ++i) {
    walker1.doSweep();
    walker2.doSweep();

    EXPECT_NEAR(walker1.getAcceptanceProbability(), walker2.getAcceptanceProbability(), 5e-7);
    EXPECT_NEAR(walker1.get_MC_log_weight(), walker2.get_MC_log_weight(), 5e-7);
    EXPECT_EQ(walker1.get_sign(), walker2.get_sign());
  }
}
