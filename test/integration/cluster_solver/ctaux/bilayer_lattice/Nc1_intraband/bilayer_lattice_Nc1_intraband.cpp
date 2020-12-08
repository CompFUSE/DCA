// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// No-change test for CT-AUX.
// Bilayer lattice with only intraband interaction.

#include <complex>
#include <iostream>
#include <string>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/function/util/difference.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_loop/dca_loop_data.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/ctaux_cluster_solver.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/models/analytic_hamiltonians/bilayer_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/testing/dca_mpi_test_environment.hpp"
#include "dca/testing/minimalist_printer.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"

constexpr bool update_baseline = false;
dca::testing::DcaMpiTestEnvironment* dca_test_env;

TEST(bilayerLattice_Nc1_intraband, Self_Energy) {
  using RngType = dca::math::random::StdRandomWrapper<std::mt19937_64>;
  using DcaPointGroupType = dca::phys::domains::D4;
  using LatticeType = dca::phys::models::bilayer_lattice<DcaPointGroupType>;
  using ModelType = dca::phys::models::TightBindingModel<LatticeType>;
  using Threading = dca::parallel::NoThreading;
  using ParametersType =
      dca::phys::params::Parameters<dca::testing::DcaMpiTestEnvironment::ConcurrencyType, Threading,
                                    dca::profiling::NullProfiler, ModelType, RngType,
                                    dca::phys::solver::CT_AUX>;
  using DcaDataType = dca::phys::DcaData<ParametersType>;
  using QmcSolverType =
      dca::phys::solver::CtauxClusterSolver<dca::linalg::CPU, ParametersType, DcaDataType>;

  using w = dca::func::dmn_0<dca::phys::domains::frequency_domain>;
  using b = dca::func::dmn_0<dca::phys::domains::electron_band_domain>;
  using s = dca::func::dmn_0<dca::phys::domains::electron_spin_domain>;
  using nu = dca::func::dmn_variadic<b, s>;  // orbital-spin index
  using k_DCA = dca::func::dmn_0<dca::phys::domains::cluster_domain<
      double, LatticeType::DIMENSION, dca::phys::domains::CLUSTER,
      dca::phys::domains::MOMENTUM_SPACE, dca::phys::domains::BRILLOUIN_ZONE>>;

  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
    dca::util::GitVersion::print();
    dca::util::Modules::print();
  }

  ParametersType parameters(dca::util::GitVersion::string(), dca_test_env->concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(dca_test_env->input_file_name);
  parameters.update_model();
  parameters.update_domains();

  dca::phys::DcaLoopData<ParametersType> dca_loop_data;

  DcaDataType dca_data_imag(parameters);
  dca_data_imag.initialize();

  // Read and broadcast ED data.
  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
    dca::io::HDF5Reader reader;
    reader.open_file(DCA_SOURCE_DIR
                     "/test/integration/cluster_solver/ctaux/bilayer_lattice/"
                     "Nc1_intra_plus_interband/data.ED.hdf5");
    reader.open_group("functions");
    // reader.execute(dca_data_imag.Sigma);
    reader.execute(dca_data_imag.G0_k_w_cluster_excluded);
    reader.execute(dca_data_imag.G0_r_t_cluster_excluded);
    reader.close_file();
  }

  // dca_test_env->concurrency.broadcast(dca_data_imag.Sigma);
  dca_test_env->concurrency.broadcast(dca_data_imag.G0_k_w_cluster_excluded);
  dca_test_env->concurrency.broadcast(dca_data_imag.G0_r_t_cluster_excluded);

  // dca::func::function<std::complex<double>, dca::func::dmn_variadic<nu, nu, k_DCA, w> >
  //   Sigma_ED(dca_data_imag.Sigma);

  // Do one QMC iteration
  QmcSolverType qmc_solver(parameters, dca_data_imag);
  qmc_solver.initialize(0);
  qmc_solver.integrate();
  qmc_solver.finalize(dca_loop_data);

  dca::func::function<std::complex<double>, dca::func::dmn_variadic<nu, nu, k_DCA, w>> Sigma_QMC(
      dca_data_imag.Sigma);

  // Read QMC self-energy from check_data file and compare it with the newly
  // computed QMC self-energy.
  const std::string filename = DCA_SOURCE_DIR
      "/test/integration/cluster_solver/ctaux/bilayer_lattice/Nc1_intraband/data.ED.hdf5";
  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
    if (!update_baseline) {
      dca::func::function<std::complex<double>, dca::func::dmn_variadic<nu, nu, k_DCA, w>> Sigma_QMC_check(
          "Self_Energy");
      dca::io::HDF5Reader reader;
      reader.open_file(filename);
      reader.open_group("functions");
      reader.execute(Sigma_QMC_check);
      reader.close_file();

      auto diff = dca::func::util::difference(Sigma_QMC_check, Sigma_QMC);
      EXPECT_GT(1e-10, diff.l_inf);
    }
    else {
      // Write results
      std::cout << "\nProcessor " << dca_test_env->concurrency.id() << " is writing data "
                << std::endl;

      dca::io::HDF5Writer writer;
      writer.open_file(filename);
      writer.open_group("functions");
      Sigma_QMC.set_name("Self_Energy");
      writer.execute(Sigma_QMC);
      writer.close_group();
      writer.close_file();
    }
    std::cout << "\nDCA main ending.\n" << std::endl;
  }
}

int main(int argc, char** argv) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  dca::parallel::MPIConcurrency concurrency(argc, argv);
  dca_test_env = new dca::testing::DcaMpiTestEnvironment(
      concurrency,
      DCA_SOURCE_DIR
      "/test/integration/cluster_solver/ctaux/bilayer_lattice/Nc1_intraband/"
      "input.bilayer_lattice_Nc1_intraband.json");
  ::testing::AddGlobalTestEnvironment(dca_test_env);

  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

  if (dca_test_env->concurrency.id() != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();

  return result;
}
