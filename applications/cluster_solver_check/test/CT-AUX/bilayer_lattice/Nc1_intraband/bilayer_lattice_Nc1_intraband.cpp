// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// No-change test for CT-AUX.
// Bilayer lattice with only intraband interaction.

#define DCA_WITH_REDUCED_VERTEX_FUNCTION

#include <complex>
#include <iostream>
#include <string>

#include "gtest/gtest.h"

#include "dca/parallel/pthreading/pthreading.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/models/analytic_hamiltonians/bilayer_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/testing/dca_mpi_test_environment.hpp"
#include "dca/testing/minimalist_printer.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"
#include "comp_library/function_library/include_function_library.h"
#include "comp_library/IO_library/HDF5/HDF5.hpp"
#include "comp_library/IO_library/JSON/JSON.hpp"
#include "comp_library/profiler_library/profilers/null_profiler.hpp"
#include "phys_library/DCA+_data/DCA_data.h"
#include "phys_library/DCA+_loop/DCA_loop_data.hpp"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_cluster_solver.h"
#include "phys_library/domains/cluster/symmetries/point_groups/2D/2D_square.h"
#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain.h"

dca::testing::DcaMpiTestEnvironment* dca_test_env;

using namespace DCA;

TEST(bilayerLattice_Nc1_intraband, Self_Energy) {
  using RngType = dca::math::random::StdRandomWrapper<std::ranlux48_base>;
  using DcaPointGroupType = D4;
  using LatticeType = dca::phys::models::bilayer_lattice<DcaPointGroupType>;
  using ModelType = dca::phys::models::TightBindingModel<LatticeType>;
  using Threading = dca::parallel::Pthreading;
  using ParametersType =
      dca::phys::params::Parameters<dca::testing::DcaMpiTestEnvironment::ConcurrencyType, Threading,
                                    PROFILER::NullProfiler, ModelType, RngType, CT_AUX_CLUSTER_SOLVER>;
  using DcaDataType = DCA_data<ParametersType>;
  using QmcSolverType =
      cluster_solver<CT_AUX_CLUSTER_SOLVER, LIN_ALG::CPU, ParametersType, DcaDataType>;

  using w = dmn_0<frequency_domain>;
  using b = dmn_0<electron_band_domain>;
  using s = dmn_0<electron_spin_domain>;
  using nu = dmn_variadic<b, s>;  // orbital-spin index
  using k_DCA =
      dmn_0<cluster_domain<double, LatticeType::DIMENSION, CLUSTER, MOMENTUM_SPACE, BRILLOUIN_ZONE>>;

  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
    dca::util::GitVersion::print();
    dca::util::Modules::print();
  }

  ParametersType parameters(dca::util::GitVersion::string(), dca_test_env->concurrency);
  parameters.read_input_and_broadcast<IO::reader<IO::JSON>>(dca_test_env->input_file_name);
  parameters.update_model();
  parameters.update_domains();

  DCA_loop_data<ParametersType> dca_loop_data;

  DcaDataType dca_data_imag(parameters);
  dca_data_imag.initialize();

  // Read and broadcast<IO::JSON> ED data
  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
    IO::reader<IO::HDF5> reader;
    reader.open_file(DCA_SOURCE_DIR
                     "/applications/cluster_solver_check/test/CT-AUX/bilayer_lattice/Nc1_interband/"
                     "data.ED.hdf5");
    reader.open_group("functions");
    // reader.execute(dca_data_imag.Sigma);
    reader.execute(dca_data_imag.G0_k_w_cluster_excluded);
    reader.execute(dca_data_imag.G0_r_t_cluster_excluded);
    reader.close_file();
  }

  // dca_test_env->concurrency.broadcast(dca_data_imag.Sigma);
  dca_test_env->concurrency.broadcast(dca_data_imag.G0_k_w_cluster_excluded);
  dca_test_env->concurrency.broadcast(dca_data_imag.G0_r_t_cluster_excluded);

  // FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_DCA, w> >
  //   Sigma_ED(dca_data_imag.Sigma);

  // Do one QMC iteration
  QmcSolverType qmc_solver(parameters, dca_data_imag);
  qmc_solver.initialize(1);
  qmc_solver.integrate();
  qmc_solver.finalize(dca_loop_data);

  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_DCA, w>> Sigma_QMC(dca_data_imag.Sigma);

  // Read QMC self-energy from check_data file and compare it with the newly
  // computed QMC self-energy.
  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_DCA, w>> Sigma_QMC_check(
        "Self_Energy");
    IO::reader<IO::HDF5> reader;
    reader.open_file(DCA_SOURCE_DIR
                     "/applications/cluster_solver_check/test/CT-AUX/bilayer_lattice/Nc1_interband/"
                     "check_data.QMC.hdf5");
    reader.open_group("functions");
    reader.execute(Sigma_QMC_check);
    reader.close_file();

    for (int w_ind = 0; w_ind < w::dmn_size(); ++w_ind) {
      for (int k_ind = 0; k_ind < k_DCA::dmn_size(); ++k_ind) {
        for (int nu_ind_2 = 0; nu_ind_2 < nu::dmn_size(); ++nu_ind_2) {
          for (int nu_ind_1 = 0; nu_ind_1 < nu::dmn_size(); ++nu_ind_1) {
            EXPECT_NEAR(Sigma_QMC_check(nu_ind_1, nu_ind_2, k_ind, w_ind).real(),
                        Sigma_QMC(nu_ind_1, nu_ind_2, k_ind, w_ind).real(), 1.e-12);
            EXPECT_NEAR(Sigma_QMC_check(nu_ind_1, nu_ind_2, k_ind, w_ind).imag(),
                        Sigma_QMC(nu_ind_1, nu_ind_2, k_ind, w_ind).imag(), 1.e-12);
          }
        }
      }
    }
  }

  // Write results
  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.last()) {
    std::cout << "\nProcessor " << dca_test_env->concurrency.id() << " is writing data " << std::endl;
    IO::writer<IO::HDF5> writer;
    writer.open_file("output.hdf5");
    writer.open_group("functions");
    Sigma_QMC.get_name() = "Self_Energy";
    writer.execute(Sigma_QMC);
    writer.close_group();
    writer.close_file();
    std::cout << "\nDCA main ending.\n" << std::endl;
  }
}
int main(int argc, char** argv) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  dca_test_env =
      new dca::testing::DcaMpiTestEnvironment(argc, argv, DCA_SOURCE_DIR
                                              "/applications/cluster_solver_check/test/"
                                              "CT-AUX/bilayer_lattice/Nc1_interband/"
                                              "input.bilayer_lattice_Nc1_interband.json");
  ::testing::AddGlobalTestEnvironment(dca_test_env);

  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

  if (dca_test_env->concurrency.id() != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();

  return result;
}
