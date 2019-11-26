// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// No-change test for a DCA+ calculation using the threaded CT-AUX cluster solver.
// It runs a simulation of a tight-binding model on 2D square lattice.

#include <iostream>
#include <string>

#include "gtest/gtest.h"

#include "dca/config/cmake_options.hpp"
#include "dca/config/threading.hpp"
//
#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/function/util/difference.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_loop/dca_loop.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/ctaux_cluster_solver.hpp"
#include "dca/phys/dca_step/cluster_solver/stdthread_qmci/stdthread_qmci_cluster_solver.hpp"

#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/testing/minimalist_printer.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"

constexpr bool update_baseline = false;

TEST(dca_sp_DCAplus_thread, Self_energy) {
#ifdef ATTACH_DEBUG
  std::cout << "Please press <return> after attaching debugger" << std::endl;
  char c;
  std::cin >> c;
#endif  // ATTACH_DEBUG

  using RngType = dca::math::random::StdRandomWrapper<std::mt19937_64>;
  using DcaPointGroupType = dca::phys::domains::D4;
  using LatticeType = dca::phys::models::square_lattice<DcaPointGroupType>;
  using ModelType = dca::phys::models::TightBindingModel<LatticeType>;
  using Concurrency = dca::parallel::NoConcurrency;
  using ParametersType =
      dca::phys::params::Parameters<Concurrency, Threading, dca::profiling::NullProfiler,
                                    ModelType, RngType, dca::phys::solver::CT_AUX>;
  using DcaDataType = dca::phys::DcaData<ParametersType>;
  using ClusterSolverBaseType =
      dca::phys::solver::CtauxClusterSolver<dca::linalg::CPU, ParametersType, DcaDataType>;

  using ClusterSolverType = dca::phys::solver::StdThreadQmciClusterSolver<ClusterSolverBaseType>;

  using DcaLoopType = dca::phys::DcaLoop<ParametersType, DcaDataType, ClusterSolverType, dca::DistType::NONE>;

  using w = dca::func::dmn_0<dca::phys::domains::frequency_domain>;
  using b = dca::func::dmn_0<dca::phys::domains::electron_band_domain>;
  using s = dca::func::dmn_0<dca::phys::domains::electron_spin_domain>;
  using nu = dca::func::dmn_variadic<b, s>;  // orbital-spin index
  using k_DCA = dca::func::dmn_0<dca::phys::domains::cluster_domain<
      double, LatticeType::DIMENSION, dca::phys::domains::CLUSTER,
      dca::phys::domains::MOMENTUM_SPACE, dca::phys::domains::BRILLOUIN_ZONE>>;

  Concurrency concurrency(0, nullptr);

  dca::util::GitVersion::print();
  dca::util::Modules::print();
  dca::config::CMakeOptions::print();

  std::cout << "\n"
            << "********************************************************************************\n"
            << "**********                     DCA(+) Calculation                     **********\n"
            << "********************************************************************************\n"
            << "\n"
            << "Start time : " << dca::util::print_time() << "\n"
            << "\n"
            << std::endl;

  ParametersType parameters(dca::util::GitVersion::string(), concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(
      DCA_SOURCE_DIR "/test/system-level/dca/input.dca_sp_DCA+_thread_test.json");
  parameters.update_model();
  parameters.update_domains();

  DcaDataType dca_data(parameters);
  dca_data.initialize();

  DcaLoopType dca_loop(parameters, dca_data, concurrency);
  dca_loop.initialize();
  dca_loop.execute();
  dca_loop.finalize();

  std::cout << "\nChecking data.\n" << std::endl;
  auto& Sigma_DCA = dca_data.Sigma;

  // Read QMC self-energy from check_data file and compare it with the newly
  // computed QMC self-energy.
  const std::string filename =
      DCA_SOURCE_DIR "/test/system-level/dca/check_data.dca_sp_DCA+_thread_test.hdf5";

  if (!update_baseline) {
    dca::func::function<std::complex<double>, dca::func::dmn_variadic<nu, nu, k_DCA, w>> Sigma_QMC_check(
        "Self_Energy");
    dca::io::HDF5Reader reader;
    reader.open_file(filename);
    reader.open_group("functions");
    reader.execute(Sigma_QMC_check);
    reader.close_file();

    auto diff = dca::func::util::difference(Sigma_QMC_check, Sigma_DCA);
    EXPECT_GT(1e-10, diff.l_inf);
  }
  else {
    // Write results
    std::cout << "\nWriting data " << std::endl;

    dca::io::HDF5Writer writer;
    writer.open_file(filename);
    writer.open_group("functions");
    Sigma_DCA.set_name("Self_Energy");
    writer.execute(Sigma_DCA);
    writer.close_group();
    writer.close_file();
  }
  std::cout << "\nDCA main ending.\n" << std::endl;
}
