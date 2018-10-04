// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// No-change test for the stdthread solver wrapper. The base solver is CT-AUX and the model is
// a square lattice with nearest neighbours hopping. Two and single particles Green's function are
// tested.

#include <iostream>
#include <string>

#include "gtest/gtest.h"

#include "dca/function/function.hpp"
#include "dca/function/util/difference.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_loop/dca_loop_data.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/ctaux_cluster_solver.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_cluster_solver.hpp"
#include "dca/phys/domains/cluster/momentum_exchange_domain.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_exchange_domain.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/stdthread/stdthread.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"

constexpr bool update_baseline = false;

const std::string input_dir = DCA_SOURCE_DIR "/test/integration/cluster_solver/stdthread_qmci/";

using Concurrency = dca::parallel::NoConcurrency;
using RngType = dca::math::random::StdRandomWrapper<std::mt19937_64>;
using Lattice = dca::phys::models::square_lattice<dca::phys::domains::D4>;
using Model = dca::phys::models::TightBindingModel<Lattice>;
using Threading = dca::parallel::stdthread;
using Parameters = dca::phys::params::Parameters<Concurrency, Threading, dca::profiling::NullProfiler,
                                                 Model, RngType, dca::phys::solver::CT_AUX>;
using Data = dca::phys::DcaData<Parameters>;
using BaseSolver = dca::phys::solver::CtauxClusterSolver<dca::linalg::CPU, Parameters, Data>;
using QmcSolver = dca::phys::solver::StdThreadQmciClusterSolver<BaseSolver>;

void performTest(const std::string& input, const std::string& baseline) {
  Concurrency concurrency(0, nullptr);
  if (concurrency.id() == concurrency.first()) {
    dca::util::GitVersion::print();
    dca::util::Modules::print();
  }

  Parameters parameters(dca::util::GitVersion::string(), concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(input_dir + input);

  static bool update_model = true;
  if (update_model) {
    parameters.update_model();
    parameters.update_domains();
  }
  else {
    // The computed frequency and momentum exchange is stored in these two domains, which would not
    // be updated otherwise.
    dca::phys::domains::FrequencyExchangeDomain::initialize(parameters);
    dca::phys::domains::MomentumExchangeDomain::initialize(parameters);
  }
  update_model = false;

  // Initialize data with G0 computation.
  Data data(parameters);
  data.initialize();

  // Do one integration step.
  QmcSolver qmc_solver(parameters, data);
  qmc_solver.initialize(0);
  qmc_solver.integrate();
  dca::phys::DcaLoopData<Parameters> loop_data;
  qmc_solver.finalize(loop_data);

  if (!update_baseline) {
    // Read and confront with previous run.
    if (concurrency.id() == 0) {
      Data::SpGreensFunction G_k_w_check(data.G_k_w.get_name());
      Data::TpGreensFunction G4_check(data.get_G4().get_name());
      G_k_w_check.set_name(data.G_k_w.get_name());
      dca::io::HDF5Reader reader;
      reader.open_file(input_dir + baseline);
      reader.open_group("functions");
      reader.execute(G_k_w_check);
      reader.execute(G4_check);
      reader.close_group(), reader.close_file();

      const auto err_g = dca::func::util::difference(G_k_w_check, data.G_k_w);
      const auto err_g4 = dca::func::util::difference(G4_check, data.get_G4());

      EXPECT_GE(5e-7, err_g.l_inf);
      EXPECT_GE(5e-7, err_g4.l_inf);
    }
  }
  else {
    //  Write results
    if (concurrency.id() == concurrency.first()) {
      dca::io::HDF5Writer writer;
      writer.open_file(input_dir + baseline);
      writer.open_group("functions");
      writer.execute(data.G_k_w);
      writer.execute(data.get_G4());
      writer.close_group(), writer.close_file();
    }
  }
}

// Test with walk and accumulation running on the same thread.
TEST(StdhreadCtauxTest, Shared) {
  performTest("stdthread_ctaux_tp_test_shared_input.json",
              "stdthread_ctaux_tp_test_shared_baseline.hdf5");
}

// Test with walk and accumulation running on different threads.
TEST(StdhreadCtauxTest, NonShared) {
  performTest("stdthread_ctaux_tp_test_nonshared_input.json",
              "stdthread_ctaux_tp_test_nonshared_baseline.hdf5");
}
