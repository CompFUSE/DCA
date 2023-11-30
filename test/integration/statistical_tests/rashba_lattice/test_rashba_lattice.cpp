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
#include <type_traits>
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
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/g0_interpolation/g0_interpolation_gpu.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/models/analytic_hamiltonians/rashba_hubbard.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/parameters.hpp"
#include "dca/profiling/null_profiler.hpp"
#include "dca/testing/minimalist_printer.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"
#include "rashba_lattice_setup.hpp"
#include "dca/math/statistical_testing/function_cut.hpp"
#include "dca/math/statistical_testing/statistical_testing.hpp"
#include "dca/testing/dca_mpi_test_environment.hpp"

constexpr bool update_baseline = false;

#ifdef DCA_HAVE_ADIOS2
#include "dca/io/adios2/adios2_reader.hpp"
#include "dca/io/adios2/adios2_writer.hpp"
adios2::ADIOS* adios2_ptr;
#endif

dca::testing::DcaMpiTestEnvironment* dca_test_env;

// #ifdef DCA_HAVE_MPI
// #include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"
// dca::parallel::MPIConcurrency* concurrency_ptr;
// #else
// #include "dca/parallel/no_concurrency/no_concurrency.hpp"
// dca::parallel::NoConcurrency* concurrency_ptr;
// #endif

using dca::linalg::DeviceType;

struct RashbaLatticeIntegrationTest : public ::testing::Test {
  using SCALAR = std::complex<double>;
  template <DeviceType DEVICE>
  using IntegrationSetupBare =
      dca::testing::IntegrationSetupBare<SCALAR, DEVICE, dca::testing::DcaMpiTestEnvironment::ConcurrencyType,
                                         dca::testing::LatticeRashba, dca::ClusterSolverId::CT_AUX>;
  virtual void reallySetUp(Concurrency* concurrency) {
    cpu_setup = std::make_unique<IntegrationSetupBare<dca::linalg::CPU>>(concurrency);
    gpu_setup = std::make_unique<IntegrationSetupBare<dca::linalg::GPU>>(concurrency);
  }

  virtual void TearDown() {
    cpu_setup.reset();
    gpu_setup.reset();
  }

  using CPUSetup = IntegrationSetupBare<dca::linalg::CPU>;
  using GPUSetup = IntegrationSetupBare<dca::linalg::GPU>;

  std::unique_ptr<CPUSetup> cpu_setup;
  std::unique_ptr<GPUSetup> gpu_setup;
};

TEST_F(RashbaLatticeIntegrationTest, SelfEnergy) {
  auto& concurrency = dca_test_env->concurrency;
  if (concurrency.id() == concurrency.first()) {
    dca::util::GitVersion::print();
    dca::util::Modules::print();
  }

  this->reallySetUp(&concurrency);

  // Do one QMC iteration
  using QMCSolverGPU = typename GPUSetup::ThreadedSolver;

  QMCSolverGPU qmc_solver_gpu(this->gpu_setup->parameters_, *(this->gpu_setup->data_), nullptr);
  qmc_solver_gpu.initialize(0);

  using QMCSolverCPU = typename RashbaLatticeIntegrationTest::CPUSetup::ThreadedSolver;
  QMCSolverCPU qmc_solver_cpu(this->cpu_setup->parameters_, *(this->cpu_setup->data_), nullptr);
  qmc_solver_cpu.initialize(0);

  using NuDmn = RashbaLatticeIntegrationTest::GPUSetup::NuDmn;
  using KDmnDCA = RashbaLatticeIntegrationTest::GPUSetup::KDmnDCA;
  using WDmn = RashbaLatticeIntegrationTest::GPUSetup::WDmn;

  auto& cpu_parameters = this->cpu_setup->parameters_;
  auto& gpu_parameters = this->gpu_setup->parameters_;

  // compare what?
  dca::linalg::Matrix<Scalar, dca::linalg::CPU> gpu_akima_coefficients_CPU;
  // G0Interpolation<dca::linalg::GPU, Parameters> g0_gpu(0, gpu_parameters);
  // gpu_akima_coefficients_CPU.setAsyncw(qmc_solver_gpu.getG0().getAkimaCoefficients(), 0, 0);
  dca::linalg::Matrix<Scalar, dca::linalg::CPU> cpu_akima_coefficients_CPU;
  auto& gpu_akima_coeff_func = qmc_solver_gpu.getG0().getAkimaCoefficients();
  auto& cpu_akima_coeff_func = qmc_solver_cpu.getG0().getAkimaCoefficients();
  EXPECT_EQ(cpu_akima_coeff_func, gpu_akima_coeff_func);
  dca::linalg::Vector<Scalar, dca::linalg::CPU> gpu_G0_coeff_on_dev;
  gpu_G0_coeff_on_dev.setAsync(qmc_solver_gpu.getG0().getG0Dev(), 0);
  dca::linalg::Vector<Scalar, dca::linalg::CPU> host_coeff;
  std::copy_n(cpu_akima_coeff_func.data(), host_coeff.size(), host_coeff.data());
  dca::linalg::util::syncStream(0, 0);
  EXPECT_EQ(gpu_G0_coeff_on_dev, host_coeff);

  auto& gpu_gkw = this->gpu_setup->data_->G_k_w;
  auto& cpu_gkw = this->cpu_setup->data_->G_k_w;
  auto diff_Gkw = dca::func::util::difference(gpu_gkw, cpu_gkw);

  dca::phys::DcaLoopData<decltype(this->gpu_setup->parameters_)> dca_loop_data_gpu;
  qmc_solver_gpu.integrate();
  auto L2_sigma_diff_gpu = qmc_solver_gpu.finalize(dca_loop_data_gpu);

  dca::phys::DcaLoopData<decltype(this->cpu_setup->parameters_)> dca_loop_data_cpu;
  qmc_solver_cpu.integrate();
  auto L2_sigma_diff_cpu = qmc_solver_cpu.finalize(dca_loop_data_cpu);

  std::cout << "L2s: " << L2_sigma_diff_gpu << " " << L2_sigma_diff_cpu << '\n';

  this->gpu_setup->performLatticeMapping();
  this->cpu_setup->performLatticeMapping();

  dca::func::function<std::complex<double>, dca::func::dmn_variadic<NuDmn, NuDmn, KDmnDCA, WDmn>>
      Sigma_QMC_gpu(this->gpu_setup->data_->Sigma);

  dca::func::function<std::complex<double>, dca::func::dmn_variadic<NuDmn, NuDmn, KDmnDCA, WDmn>>
      Sigma_QMC_cpu(this->cpu_setup->data_->Sigma);
  dca::linalg::util::syncStream(0, 0);

  auto diff = dca::func::util::difference(Sigma_QMC_gpu, Sigma_QMC_cpu);
  EXPECT_LT(diff.l_inf, 1E-10);

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
  // #ifdef DCA_HAVE_MPI
  //   dca::parallel::MPIConcurrency concurrency(argc, argv);
  //   concurrency_ptr = &concurrency;
  // #else
  //   dca::parallel::NoConcurrency concurrency(argc, argv);
  //   concurrency_ptr = &concurrency;
  // #endif
  ::testing::InitGoogleTest(&argc, argv);

  dca::parallel::MPIConcurrency concurrency(argc, argv);
  dca_test_env =
      new dca::testing::DcaMpiTestEnvironment(concurrency, "rashba_lattice_stat_test.json");
  testing::AddGlobalTestEnvironment(dca_test_env);

  dca::linalg::util::initializeMagma();

#ifdef DCA_HAVE_ADIOS2
  // ADIOS expects MPI_COMM pointer or nullptr
  // std::string dummy("");
  // adios2::ADIOS adios(dummy, concurrency.get(), false);
  // adios2_ptr = &adios;
#endif

  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  // delete listeners.Release(listeners.default_result_printer());
  // listeners.Append(new dca::testing::MinimalistPrinter);
  if (dca_test_env->concurrency.id() != 0) {
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new dca::testing::MinimalistPrinter);
  }

  int result = RUN_ALL_TESTS();
  return result;
}
