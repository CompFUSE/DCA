// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter W. Doak (doakpw@ornl.gov)
//
// Comparison test for complex Rashba-Hubbardd model with CT-AUX.

#include "dca/config/profiler.hpp"
#include <complex>
#include <iostream>
#include <string>
#include <dca/function/util/difference.hpp>

#include "dca/testing/gtest_h_w_warning_blocking.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/io/hdf5/hdf5_reader.hpp"
#include "dca/io/hdf5/hdf5_writer.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/parallel/no_threading/no_threading.hpp"

using Scalar = std::complex<double>;

#include "test/mock_mcconfig.hpp"
namespace dca {
namespace config {
using McOptions = MockMcOptions<Scalar>;
}  // namespace config
}  // namespace dca

#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_loop/dca_loop_data.hpp"
#include "dca/config/profiler.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/ctaux_cluster_solver.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/models/analytic_hamiltonians/rashba_hubbard.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/testing/dca_mpi_test_environment.hpp"
#include "dca/testing/minimalist_printer.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"
#include "rashba_lattice_setup.hpp"
#include "dca/math/statistical_testing/function_cut.hpp"
#include "dca/math/statistical_testing/statistical_testing.hpp"

constexpr bool update_baseline = false;

#ifdef DCA_HAVE_ADIOS2
adios2::ADIOS* adios_ptr;
#endif

#ifdef DCA_HAVE_MPI
#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"
dca::parallel::MPIConcurrency* concurrency_ptr;
#else
#include "dca/parallel/no_concurrency/no_concurrency.hpp"
dca::parallel::NoConcurrency* concurrency_ptr;
#endif

using dca::linalg::DeviceType;

template <typename SCALAR>
struct RashbaLatticeIntegrationTest : public ::testing::Test {
  template <DeviceType DEVICE>

  using IntegrationSetupBare =
      dca::testing::IntegrationSetupBare<SCALAR, DEVICE, dca::testing::LatticeRashba,
                                         dca::ClusterSolverId::CT_AUX>;
  virtual void SetUp() {
    host_setup.SetUp();
    gpu_setup.SetUp();
  }

  virtual void TearDown() {}
  IntegrationSetupBare<dca::linalg::CPU> host_setup;
  IntegrationSetupBare<dca::linalg::GPU> gpu_setup;
};

using TestTypes = ::testing::Types<std::complex<double>>;
TYPED_TEST_CASE(RashbaLatticeIntegrationTest, TestTypes);

TYPED_TEST(RashbaLatticeIntegrationTest, SelfEnergy) {
  if (concurrency_ptr->id() == concurrency_ptr->first()) {
    dca::util::GitVersion::print();
    dca::util::Modules::print();
  }

  dca::phys::DcaLoopData<decltype(this->gpu_setup.parameters_)> dca_loop_data_gpu;

  // dca::func::function<std::complex<double>, dca::func::dmn_variadic<nu, nu, k_DCA, w> >
  //   Sigma_ED(dca_data_imag.Sigma);

  // Do one QMC iteration
  using QMCSolverGPU = typename decltype(this->gpu_setup)::ThreadedSolver;
  QMCSolverGPU qmc_solver_gpu(this->gpu_setup.parameters_, *(this->gpu_setup.data_), nullptr);
  qmc_solver_gpu.initialize(0);
  qmc_solver_gpu.integrate();
  qmc_solver_gpu.finalize(dca_loop_data_gpu);

  using NuDmn = typename decltype(this->gpu_setup)::NuDmn;
  using KDmnDCA = typename decltype(this->gpu_setup)::KDmnDCA;
  using WDmn = typename decltype(this->gpu_setup)::WDmn;

  dca::func::function<std::complex<double>, dca::func::dmn_variadic<NuDmn, NuDmn, KDmnDCA, WDmn>> Sigma_QMC(
      this->gpu_setup.data_->Sigma);

  // // Read QMC self-energy from check_data file and compare it with the newly
  // // computed QMC self-energy.
  // const std::string filename = DCA_SOURCE_DIR
  //     "/test/integration/cluster_solver/ctaux/rashba_lattice/rashba_check_data.QMC.hdf5";
  // if (dca_test_env->concurrency.id() == dca_test_env->concurrency.first()) {
  //   if (!update_baseline) {
  //     dca::func::function<std::complex<double>, dca::func::dmn_variadic<nu, nu, k_DCA, w>> Sigma_QMC_check(
  //         "SelfEnergy");
  //     dca::io::HDF5Reader reader;
  //     reader.open_file(filename);
  //     reader.open_group("functions");
  //     reader.execute(Sigma_QMC_check);
  //     reader.close_file();

  //     auto diff = dca::func::util::difference(Sigma_QMC_check, Sigma_QMC);
  //     EXPECT_GT(1e-10, diff.l_inf);
  //   }
  //   else {
  //     // Write results
  //     std::cout << "\nProcessor " << dca_test_env->concurrency.id() << " is writing data "
  //               << std::endl;

  //     dca::io::HDF5Writer writer;
  //     writer.open_file(filename);
  //     writer.open_group("functions");
  //     Sigma_QMC.set_name("Self_Energy");
  //     writer.execute(Sigma_QMC);
  //     writer.close_group();
  //     writer.close_file();
  //   }
  //   std::cout << "\nDCA main ending.\n" << std::endl;
  // }
}

int main(int argc, char** argv) {
#ifdef DCA_HAVE_MPI
  dca::parallel::MPIConcurrency concurrency(argc, argv);
  concurrency_ptr = &concurrency;
#else
  dca::parallel::NoConcurrency concurrency(argc, argv);
  concurrency_ptr = &concurrency;
#endif

  dca::linalg::util::initializeMagma();

#ifdef DCA_HAVE_ADIOS2
  // ADIOS expects MPI_COMM pointer or nullptr
  adios2::ADIOS adios("", concurrency_ptr->get(), false);
  adios_ptr = &adios;
#endif
  ::testing::InitGoogleTest(&argc, argv);

  // ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  // delete listeners.Release(listeners.default_result_printer());
  // listeners.Append(new dca::testing::MinimalistPrinter);

  int result = RUN_ALL_TESTS();
  return result;
}
