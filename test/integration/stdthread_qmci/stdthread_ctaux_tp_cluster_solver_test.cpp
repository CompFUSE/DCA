// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// No-change test for MC posix wrapper.

#include <iostream>
#include <string>

#include "gtest/gtest.h"

#include "dca/function/function.hpp"
#include "dca/function/function_utils.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_loop/dca_loop_data.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/ctaux_cluster_solver.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_cluster_solver.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/pthreading/pthreading.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/testing/minimalist_printer.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"

constexpr bool UPDATE_RESULTS = false;

const std::string input_dir = DCA_SOURCE_DIR "/test/integration/stdthread_qmci/";

using Concurrency = dca::parallel::NoConcurrency;
using RngType = dca::math::random::StdRandomWrapper<std::mt19937_64>;
using Lattice = dca::phys::models::square_lattice<dca::phys::domains::D4>;
using Model = dca::phys::models::TightBindingModel<Lattice>;
using Threading = dca::parallel::Pthreading;
using Parameters = dca::phys::params::Parameters<Concurrency, Threading, dca::profiling::NullProfiler,
                                                 Model, RngType, dca::phys::solver::CT_AUX>;
using Data = dca::phys::DcaData<Parameters>;
using BaseSolver = dca::phys::solver::CtauxClusterSolver<dca::linalg::GPU, Parameters, Data>;
using QmcSolver = dca::phys::solver::StdThreadQmciClusterSolver<BaseSolver>;

TEST(PosixCtauxClusterSolverTest, G_k_w) {
  dca::linalg::util::initializeMagma();
  Concurrency concurrency(0, nullptr);
  if (concurrency.id() == concurrency.first()) {
    dca::util::GitVersion::print();
    dca::util::Modules::print();
  }

  Parameters parameters(dca::util::GitVersion::string(), concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(input_dir + "input.json");
  parameters.update_model();
  parameters.update_domains();

  // Initialize data with G0 computation.
  Data data(parameters);
  data.initialize();

  // Do one integration step.
  QmcSolver qmc_solver(parameters, data);
  qmc_solver.initialize(0);
  qmc_solver.integrate();
  dca::phys::DcaLoopData<Parameters> loop_data;
  qmc_solver.finalize(loop_data);

  if (not UPDATE_RESULTS) {
    // Read and confront with previous run.
    if (concurrency.id() == 0) {
      auto G_k_w_check = data.G_k_w;
      using DomainType = typename Data::TpGreensFunction::this_domain_type;
      dca::func::function<std::complex<double>, DomainType> G4_check(data.get_G4_k_k_w_w().get_name());
      G_k_w_check.set_name(data.G_k_w.get_name());
      dca::io::HDF5Reader reader;
      reader.open_file(input_dir + "data.hdf5");
      reader.execute(G_k_w_check);
      reader.execute(G4_check);
      reader.close_file();

      auto err_g = dca::func::utils::difference(G_k_w_check, data.G_k_w);
      auto err_g4 = dca::func::utils::difference(G4_check, data.get_G4_k_k_w_w());

      EXPECT_GE(5e-7, err_g.l_inf);
      EXPECT_GE(5e-7, err_g4.l_inf);
    }
  }
  else {
    //  Write results
    if (concurrency.id() == concurrency.first()) {
      dca::io::HDF5Writer writer;
      writer.open_file("data.hdf5");
      writer.execute(data.G_k_w);
      writer.execute(data.get_G4_k_k_w_w());
      writer.close_file();
    }
  }
}
