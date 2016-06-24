// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// Integration test for a DCA+ main loop calculation using the pthreaded CT-AUX cluster solver.
// It runs a simulation of a tight-binding model on 2D square lattice.

#include "dca/config/defines.hpp"
#ifndef DCA_HAVE_MPI
#error MPI must be supported for the dca_DCA+_mpi_test.
#endif  // DCA_HAVE_MPI

#include <iostream>
#include <string>

#include "gtest/gtest.h"

#include "dca_mpi_test_environment.hpp"
#include "minimalist_printer.hpp"
#include "dca/concurrency/parallelization_pthreads.h"
#include "dca/math_library/random_number_library/ranq2.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"
#include "comp_library/function_library/include_function_library.h"
#include "comp_library/IO_library/HDF5/HDF5.hpp"
#include "comp_library/IO_library/JSON/JSON.hpp"
#include "phys_library/DCA+_data/DCA_data.h"
#include "phys_library/DCA+_loop/DCA_loop.hpp"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_cluster_solver.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_pthread_jacket/posix_qmci_cluster_solver.h"
#include "phys_library/domains/cluster/symmetries/point_groups/2D/2D_square.h"
#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain.h"
#include "phys_library/parameters/models/analytic_hamiltonians/lattices/2D_square_lattice.h"
#include "phys_library/parameters/models/tight_binding_model.h"
#include "phys_library/parameters/Parameters.h"

dca::testing::DcaMpiTestEnvironment* dca_test_env;

namespace dca {
namespace testing {
// dca::testing::

using namespace DCA;

TEST(dca_sp_DCAplus_pthread, Self_energy) {
#ifdef ATTACH_DEBUG
  std::cout << "Please press <return> after attaching debugger" << std::endl;
  char c;
  std::cin >> c;
#endif  // ATTACH_DEBUG

  using RngType = rng::ranq2;
  using DcaPointGroupType = D4;
  using LatticeType = square_lattice<DcaPointGroupType>;
  using ModelType = tight_binding_model<LatticeType>;
  using ParametersType =
      Parameters<DcaMpiTestEnvironment::ConcurrencyType, ModelType, RngType, CT_AUX_CLUSTER_SOLVER>;
  using DcaDataType = DCA_data<ParametersType>;
  using ClusterSolverBaseType =
      cluster_solver<CT_AUX_CLUSTER_SOLVER, LIN_ALG::CPU, ParametersType, DcaDataType>;
  using ClusterSolverType = DCA::posix_qmci_integrator<ClusterSolverBaseType>;
  using DcaLoopType = DCA_loop<ParametersType, DcaDataType, ClusterSolverType>;

  using w = dmn_0<frequency_domain>;
  using b = dmn_0<electron_band_domain>;
  using s = dmn_0<electron_spin_domain>;
  using nu = dmn_variadic<b, s>;  // orbital-spin index
  using k_DCA =
      dmn_0<cluster_domain<double, LatticeType::DIMENSION, CLUSTER, MOMENTUM_SPACE, BRILLOUIN_ZONE>>;

  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
    std::cout << "\nDCA main starting.\n"
              << "MPI-world set up: " << dca_test_env->concurrency.number_of_processors()
              << " processes.\n"
              << std::endl;

    dca::util::GitVersion::print();
    dca::util::Modules::print();
  }

  ParametersType parameters(dca::util::GitVersion::string(), dca_test_env->concurrency);
  parameters.read_input_and_broadcast<IO::reader<IO::JSON>>(dca_test_env->input_file_name);
  parameters.update_model();
  parameters.update_domains();

  DcaDataType dca_data(parameters);
  dca_data.initialize();

  DcaLoopType dca_loop(parameters, dca_data, dca_test_env->concurrency);
  dca_loop.initialize();
  dca_loop.execute();
  dca_loop.finalize();

  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
    std::cout << "\nProcessor " << dca_test_env->concurrency.id() << " is checking data "
              << std::endl;

    // Read self-energy from check_data file.
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, k_DCA, w>> Sigma_check("Self_Energy");
    IO::reader<IO::HDF5> reader;
    reader.open_file(DCA_SOURCE_DIRECTORY
                     "/applications/dca/test/check_data.dca_sp_DCA+_pthread_test.hdf5");
    reader.open_group("functions");
    reader.execute(Sigma_check);
    reader.close_file();

    // Compare the computed self-energy with the expected result.
    for (int w_ind = 0; w_ind < w::dmn_size(); ++w_ind) {
      for (int k_ind = 0; k_ind < k_DCA::dmn_size(); ++k_ind) {
        for (int nu_ind_2 = 0; nu_ind_2 < nu::dmn_size(); ++nu_ind_2) {
          for (int nu_ind_1 = 0; nu_ind_1 < nu::dmn_size(); ++nu_ind_1) {
            EXPECT_NEAR(Sigma_check(nu_ind_1, nu_ind_2, k_ind, w_ind).real(),
                        dca_data.Sigma(nu_ind_1, nu_ind_2, k_ind, w_ind).real(), 1.e-12);
            EXPECT_NEAR(Sigma_check(nu_ind_1, nu_ind_2, k_ind, w_ind).imag(),
                        dca_data.Sigma(nu_ind_1, nu_ind_2, k_ind, w_ind).imag(), 1.e-12);
          }
        }
      }
    }
  }

  if (dca_test_env->concurrency.id() == dca_test_env->concurrency.last()) {
    std::cout << "\nProcessor " << dca_test_env->concurrency.id() << " is writing data " << std::endl;
    dca_loop.write();

    std::cout << "\nDCA main ending.\n" << std::endl;
  }
}

}  // testing
}  // dca

int main(int argc, char** argv) {
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  dca_test_env = new dca::testing::DcaMpiTestEnvironment(
      argc, argv,
      DCA_SOURCE_DIRECTORY "/applications/dca/test/input.dca_sp_DCA+_pthread_test.json");
  ::testing::AddGlobalTestEnvironment(dca_test_env);

  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

  if (dca_test_env->concurrency.id() != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  result = RUN_ALL_TESTS();

  return result;
}
