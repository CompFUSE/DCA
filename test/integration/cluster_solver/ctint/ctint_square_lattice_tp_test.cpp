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
// Square lattice with single band and double occupancy repulsion U.
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
#include "dca/phys/dca_step/cluster_solver/ctint/ctint_cluster_solver.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_cluster_solver.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/parallel/stdthread/stdthread.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"

const std::string input_dir =
    DCA_SOURCE_DIR "/test/integration/cluster_solver/ctint/";

constexpr bool update_baseline = false;

TEST(CtintSquareLatticeTpTest, Self_Energy) {
  using RngType = dca::math::random::StdRandomWrapper<std::ranlux48_base>;
  using Lattice = dca::phys::models::square_lattice<dca::phys::domains::D4>;
  using Model = dca::phys::models::TightBindingModel<Lattice>;
  using Concurrency = dca::parallel::NoConcurrency;
  using Threading = dca::parallel::stdthread;
  using Parameters =
      dca::phys::params::Parameters<Concurrency, Threading, dca::profiling::NullProfiler, Model,
                                    RngType, dca::phys::solver::CT_INT>;
  using Data = dca::phys::DcaData<Parameters>;
  using BaseSolverType = dca::phys::solver::CtintClusterSolver<dca::linalg::CPU, Parameters>;
  using QmcSolverType = dca::phys::solver::StdThreadQmciClusterSolver<BaseSolverType>;

  Concurrency concurrency(0, nullptr);
  dca::util::GitVersion::print();
  dca::util::Modules::print();

  Parameters parameters(dca::util::GitVersion::string(), concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(input_dir +
                                                           "square_lattice_tp_input.json");
  parameters.update_model();
  parameters.update_domains();

  // Initialize data with G0 computation.
  Data data(parameters);
  data.initialize();

  // Do one integration step.
  QmcSolverType qmc_solver(parameters, data);
  qmc_solver.initialize(0);
  qmc_solver.integrate();

  dca::phys::DcaLoopData<Parameters> loop_data;
  qmc_solver.finalize(loop_data);

  if (!update_baseline) {
    // Read and confront with previous run
    Data::TpGreensFunction G4_check(data.get_G4()[0].get_name());
    Data::SpGreensFunction G_check(data.G_k_w.get_name());
    dca::io::HDF5Reader reader;
    reader.open_file(input_dir + "square_lattice_tp_baseline.hdf5");
    reader.open_group("functions");
    reader.execute("G4", G4_check);
    reader.execute(G_check);
    reader.close_group(), reader.close_file();

    const auto diff_g = dca::func::util::difference(G_check, data.G_k_w);
    const auto diff_g4 = dca::func::util::difference(G4_check, data.get_G4()[0]);
    EXPECT_GE(5e-7, diff_g.l_inf);
    EXPECT_GE(5e-7, diff_g4.l_inf);
  }
  else {
    //  Write results
    dca::io::HDF5Writer writer;
    writer.open_file(input_dir + "square_lattice_tp_baseline.hdf5");
    writer.open_group("functions");
    writer.execute("G4", data.get_G4()[0]);
    writer.execute(data.G_k_w);
    writer.close_group(), writer.close_file();
  }
}
