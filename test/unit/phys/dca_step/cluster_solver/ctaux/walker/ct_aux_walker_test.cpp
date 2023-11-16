// Copyright (C) 2023 UT-Battelle, LLC
// Copyright (C) 2023 ETH Zurich
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter W. Doak (doakpw@ornl.gov)
//
// This file tests the CPU/GPU walkers used by the ctaux cluster solver.
// based on ct_int_walker_test.cpp

#include "dca/platform/dca_gpu.h"
#include "dca/phys/dca_step/cluster_solver/ctaux/ctaux_walker.hpp"
#include "dca/testing/gtest_h_w_warning_blocking.h"
#include "walker_wrapper.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/util/to_string.hpp"

#include "test/mock_mcconfig.hpp"
namespace dca {
namespace config {
using McOptions = MockMcOptions<std::complex<double>>;
}  // namespace config
}  // namespace dca

constexpr char rashba_input[] =
    DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/ctaux/walker/input_rashba_ctaux.json";

// template <typename Scalar>
// using CtauxWalkerTestT =
//   typename dca::testing::G0Setup<Scalar, dca::testing::LatticeRashba, dca::ClusterSolverId::CT_AUX, rashba_input>;

// Todo parameterize this test on the Lattice and get coverage of both double and std::complex<double>
template <typename SCALAR>
struct CtauxWalkerTestT : public ::testing::Test {
  using G0Setup = dca::testing::G0SetupBare<SCALAR, dca::testing::LatticeRashba,
                                            dca::ClusterSolverId::CT_AUX, rashba_input>;
  virtual void SetUp() {
    host_setup.SetUp();
    gpu_setup.SetUp();
  }

  virtual void TearDown() {}
  G0Setup host_setup;
  G0Setup gpu_setup;
};

using Scalar = std::complex<double>;
using CtauxWalkerTest = CtauxWalkerTestT<Scalar>;

using namespace dca::phys::solver;
using namespace dca::phys::solver::ctaux;
using namespace dca::addt_str_oper;

template <typename Scalar>
using Matrix = dca::linalg::Matrix<Scalar, dca::linalg::CPU>;
template <typename Scalar>
using MatrixPair = std::array<Matrix<Scalar>, 2>;

using dca::linalg::matrixop::determinantIP;
template <typename Scalar>
double computeDetRatio(MatrixPair<Scalar> a, MatrixPair<Scalar> b) {
  double res = 1;
  res *= determinantIP(a[0]) / determinantIP(b[0]);
  res *= determinantIP(a[1]) / determinantIP(b[1]);
  return res;
}

template <typename Scalar>
double determinant(MatrixPair<Scalar> a) {
  return determinantIP(a[0]) * determinantIP(a[1]);
}

// using FloatingPointTypes = ::testing::Types<std::complex<double>>;
// TYPED_TEST_CASE(CtauxWalkerTest, FloatingPointTypes);

TEST_F(CtauxWalkerTest, InsertAndRemoveVertex) {
  // using Scalar = TypeParam;

  // Setup
  std::vector<double> rng_values(1000);
  for (auto& x : rng_values)
    x = static_cast<double>(std::rand()) / RAND_MAX;
  typename CtauxWalkerTest::G0Setup::RngType cpu_rng(rng_values);
  typename CtauxWalkerTest::G0Setup::RngType gpu_rng(rng_values);

  auto& cpu_data = *this->host_setup.data_;
  auto& gpu_data = *this->gpu_setup.data_;
  using Parameters = typename CtauxWalkerTest::G0Setup::Parameters;
  auto& cpu_parameters = this->host_setup.parameters_;
  auto& gpu_parameters = this->gpu_setup.parameters_;

  G0Interpolation<dca::linalg::CPU, Parameters> g0_cpu(0, cpu_parameters);
  g0_cpu.initialize(cpu_data);
  G0Interpolation<dca::linalg::GPU, Parameters> g0_gpu(0, gpu_parameters);
  g0_gpu.initialize(gpu_data);

  using CPUWalker =
      testing::phys::solver::ctaux::CTAUXWalkerWrapper<Scalar, dca::linalg::CPU, Parameters>;
  using GPUWalker =
      testing::phys::solver::ctaux::CTAUXWalkerWrapper<Scalar, dca::linalg::GPU, Parameters>;

  // These break out the many sub calls in CtauxWalker::doStep so we can check values
  auto stepBeforeComputeGamma = [](auto& walker, auto& steps) {
    walker.get_configuration().prepare_configuration();
    walker.generate_delayed_spins(steps);
    walker.addNonInteractingSpinsToMatrices();
  };
  auto stepAfterComputeGamma = [](auto& walker, bool thermalized) {
    walker.upload_to_device();
    walker.update_N_matrix_with_Gamma_matrix();
    walker.clean_up_the_configuration();
    if (!thermalized)
      walker.getWarmUpdateExpansionOrder().addSample(
          walker.get_configuration().get_number_of_interacting_HS_spins());
  };

  CPUWalker cpu_walker(cpu_parameters, cpu_data, cpu_rng, 0);

  cpu_walker.initialize(0);
  // for (int i = 0; i < 1; i++) {
  //   cpu_walker.doSweep();
  //   cpu_walker.updateShell(i, cpu_parameters.get_warm_up_sweeps());
  // }

  // cpu_walker.markThermalized();

  int cpu_steps = 1;
  stepBeforeComputeGamma(cpu_walker, cpu_steps);

  GPUWalker gpu_walker(gpu_parameters, gpu_data, gpu_rng, 0);
  gpu_walker.initialize(0);

  // for (int i = 0; i < 1; i++) {
  //   gpu_walker.doSweep();
  //   gpu_walker.updateShell(i, gpu_parameters.get_warm_up_sweeps());
  // }

  // gpu_walker.markThermalized();

  int gpu_steps = 1;
  stepBeforeComputeGamma(gpu_walker, gpu_steps);

  // Both walkers have gone through states in doStep up to computeGamma

  // These values don't get moved so we need to do that to check them here
  dca::linalg::Matrix<Scalar, dca::linalg::CPU> G0_up_GPU;
  dca::linalg::Matrix<Scalar, dca::linalg::CPU> G0_dn_GPU;
  G0_up_GPU.setAsync(gpu_walker.getG0Up(), 0, 0);
  G0_dn_GPU.setAsync(gpu_walker.getG0Dn(), 0, 0);
  cudaStreamSynchronize(0);

  EXPECT_TRUE(G0_up_GPU == cpu_walker.getG0Up());
  EXPECT_TRUE(G0_dn_GPU == cpu_walker.getG0Dn())
      << G0_dn_GPU.toStr() << " == " << cpu_walker.getG0Dn().toStr() << '\n';

  dca::linalg::Matrix<Scalar, dca::linalg::CPU> G_up_GPU;
  dca::linalg::Matrix<Scalar, dca::linalg::CPU> G_dn_GPU;
  G_up_GPU.setAsync(gpu_walker.getGUp(), 0, 0);
  G_dn_GPU.setAsync(gpu_walker.getGDn(), 0, 0);
  cudaStreamSynchronize(0);

  EXPECT_TRUE(G_up_GPU == cpu_walker.getGUp());
  EXPECT_TRUE(G_dn_GPU == cpu_walker.getGDn())
      << G_dn_GPU.toStr() << " == " << cpu_walker.getGDn().toStr() << '\n';

  dca::linalg::Matrix<Scalar, dca::linalg::CPU> N_up_GPU;
  dca::linalg::Matrix<Scalar, dca::linalg::CPU> N_dn_GPU;
  N_up_GPU.setAsync(gpu_walker.getNUp(), 0, 0);
  N_dn_GPU.setAsync(gpu_walker.getNDn(), 0, 0);
  cudaStreamSynchronize(0);

  EXPECT_TRUE(N_up_GPU == cpu_walker.getNUp());
  EXPECT_TRUE(N_dn_GPU == cpu_walker.getNDn())
      << N_dn_GPU.toStr() << " == " << cpu_walker.getNDn().toStr() << '\n';

  // Up and down must be done separately as there is a state dependences between
  // calcExp and Read
  // start up
  cpu_walker.calcExpVandPushVertexIndices(dca::phys::e_UP);
  gpu_walker.calcExpVandPushVertexIndices(dca::phys::e_UP);

  dca::linalg::Vector<Scalar, dca::linalg::CPU> expV_CPU;
  dca::linalg::Vector<Scalar, dca::linalg::CPU> expV_GPU;
  expV_CPU.setAsync(cpu_walker.getExpV(), 0, 0);
  expV_GPU.setAsync(gpu_walker.getExpV(), 0, 0);
  cudaStreamSynchronize(0);
  EXPECT_EQ(expV_CPU, expV_GPU);

  dca::linalg::Vector<Scalar, dca::linalg::CPU> expDeltaV_CPU;
  dca::linalg::Vector<Scalar, dca::linalg::CPU> expDeltaV_GPU;
  expV_CPU.setAsync(cpu_walker.getExpDeltaV(), 0, 0);
  expV_GPU.setAsync(gpu_walker.getExpDeltaV(), 0, 0);
  cudaStreamSynchronize(0);
  EXPECT_EQ(expDeltaV_CPU, expDeltaV_GPU);

  dca::linalg::Vector<int, dca::linalg::CPU> vertex_ind_CPU;
  dca::linalg::Vector<int, dca::linalg::CPU> vertex_ind_GPU;
  vertex_ind_CPU.setAsync(cpu_walker.getVertexInd(), 0, 0);
  vertex_ind_GPU.setAsync(gpu_walker.getVertexInd(), 0, 0);
  cudaStreamSynchronize(0);
  EXPECT_EQ(vertex_ind_CPU, vertex_ind_GPU);

  cpu_walker.read_Gamma_matrices(dca::phys::e_UP);
  gpu_walker.read_Gamma_matrices(dca::phys::e_UP);
  // complete up

  // start dn
  cpu_walker.calcExpVandPushVertexIndices(dca::phys::e_DN);
  gpu_walker.calcExpVandPushVertexIndices(dca::phys::e_DN);

  expV_CPU.setAsync(cpu_walker.getExpV(), 0, 0);
  expV_GPU.setAsync(gpu_walker.getExpV(), 0, 0);
  cudaStreamSynchronize(0);

  EXPECT_EQ(expV_CPU, expV_GPU);

  expV_CPU.setAsync(cpu_walker.getExpDeltaV(), 0, 0);
  expV_GPU.setAsync(gpu_walker.getExpDeltaV(), 0, 0);
  cudaStreamSynchronize(0);
  EXPECT_EQ(expDeltaV_CPU, expDeltaV_GPU);

  vertex_ind_CPU.setAsync(cpu_walker.getVertexInd(), 0, 0);
  vertex_ind_GPU.setAsync(gpu_walker.getVertexInd(), 0, 0);
  cudaStreamSynchronize(0);
  EXPECT_EQ(vertex_ind_CPU, vertex_ind_GPU);

  cpu_walker.read_Gamma_matrices(dca::phys::e_DN);
  gpu_walker.read_Gamma_matrices(dca::phys::e_DN);

  cpu_walker.actually_download_from_device();
  gpu_walker.actually_download_from_device();

  auto& Gamma_up_CPU = cpu_walker.getGamma_up_CPU();
  auto& Gamma_dn_CPU = cpu_walker.getGamma_dn_CPU();
  // the ctaux walker moves the Gamma on and off the device into the Gamma_up/dn_CPU
  // it does a download to host as the last step before computeGamma
  auto& Gamma_up_GPU = gpu_walker.getGamma_up_CPU();
  auto& Gamma_dn_GPU = gpu_walker.getGamma_dn_CPU();

  // This difference can just not reduced the (gamma_k) / (gamma_k - 1) division
  // is just limited to this precision.
  auto diff = dca::linalg::matrixop::difference(Gamma_dn_CPU, Gamma_dn_GPU);
  EXPECT_LT(diff, 1E-15) << Gamma_dn_CPU.toStr() << " == " << Gamma_dn_GPU.toStr() << '\n';

  // EXPECT_TRUE(Gamma_up_CPU == Gamma_up_GPU);
  // EXPECT_TRUE(Gamma_dn_CPU == Gamma_dn_GPU);

  auto& cpu_configuration = cpu_walker.get_configuration();
  auto& gpu_configuration = gpu_walker.get_configuration();

  EXPECT_TRUE(cpu_configuration == gpu_configuration);

  // is the small difference between these causing the bigger difference after.
  // Gamma_dn_GPU = Gamma_dn_CPU;
  // The answer is yes
  dca::linalg::util::syncStream(0, 0);

  cpu_walker.compute_Gamma_matrices();
  gpu_walker.compute_Gamma_matrices();

  cudaStreamSynchronize(0);
  auto expected_size = std::pair<int, int>{0, 0};
  EXPECT_TRUE(Gamma_up_CPU.size() == expected_size);
  EXPECT_EQ(Gamma_up_CPU.size(), Gamma_up_GPU.size());
  EXPECT_EQ(Gamma_dn_CPU.size(), Gamma_dn_GPU.size());

  diff = dca::linalg::matrixop::difference(Gamma_dn_CPU, Gamma_dn_GPU);
  // even though both seem to be calculated on the CPU the tiny difference inthe initialized Gamma
  // becomes much larger, it it still insignificant?
  EXPECT_LT(diff, 1E-15) << Gamma_dn_CPU.toStr() << " == " << Gamma_dn_GPU.toStr() << '\n';
  stepAfterComputeGamma(cpu_walker, false);
  // cpu_walker.doStep(cpu_steps);

  stepAfterComputeGamma(gpu_walker, false);
  // gpu_walker.doStep(gpu_steps);

  // this passes only because Gamma_up is size 0 so nothing gets moved to the device.
  EXPECT_TRUE(Gamma_up_CPU == Gamma_up_GPU);

  // This should fail because Gamma_dn_CPU is moved to device and 0'd in size on CPU for Gamma_dn_GPU
  // diff = dca::linalg::matrixop::difference(Gamma_dn_CPU, Gamma_dn_GPU);
  // EXPECT_LT(diff, 1E-28);

  EXPECT_TRUE(cpu_configuration == gpu_configuration);

  auto& cpu_n_up = cpu_walker.getNUp();
  auto& gpu_n_up_dev = gpu_walker.getNUp();
  dca::linalg::Matrix<Scalar, dca::linalg::DeviceType::CPU> gpu_n_up{gpu_n_up_dev};
  EXPECT_EQ(cpu_n_up, gpu_n_up);

  dca::math::Phase<std::complex<double>> phase;

  auto process_matrix = [](auto& m) -> dca::math::Phase<std::complex<double>> {
    dca::math::Phase<std::complex<double>> phase_out;
    if (!m.nrRows())
      return phase_out;
    const auto [log_det, phase] = dca::linalg::matrixop::logDeterminant(m);
    std::cout << "phase: " << phase.getSign() << '\n';
    return phase;
  };

  auto cpu_proc_phase = process_matrix(cpu_n_up);
  auto gpu_proc_phase = process_matrix(gpu_n_up);
  EXPECT_NEAR(cpu_proc_phase.getSign().real(), gpu_proc_phase.getSign().real(), 1E-28);
  EXPECT_NEAR(cpu_proc_phase.getSign().imag(), gpu_proc_phase.getSign().imag(), 1E-28);

  auto cpu_phase = cpu_walker.get_sign();
  auto gpu_phase = gpu_walker.get_sign();
  EXPECT_NEAR(cpu_phase.real(), gpu_phase.real(), 1E-28);
  EXPECT_NEAR(cpu_phase.imag(), gpu_phase.imag(), 1E-28);

  // // *******************************
  // // Test vertex removal ***********
  // // *******************************
  // // Set rng value to select: last vertex, unused, unused, accept
  // rng.setNewValues(std::vector<Scalar>{0.95, -1, -1, 0.01});
  // MatrixPair<Scalar> old_M(walker.getM());
  // bool result = walker.tryVertexRemoval();
  // MatrixPair<Scalar> new_M(walker.getM());
  // ASSERT_EQ(true, result);
  // //  ASSERT_EQ(old_M.nrCols(), new_M.nrCols() + 2);
  // // Compute directly the new M.
  // walker.setMFromConfig();
  // MatrixPair<Scalar> direct_M(walker.getM());
  // for (int s = 0; s < 2; ++s)
  //   for (int j = 0; j < new_M[s].nrCols(); j++)
  //     for (int i = 0; i < new_M[s].nrRows(); i++)
  //       EXPECT_NEAR(direct_M[s](i, j), new_M[s](i, j), 10 * std::numeric_limits<Scalar>::epsilon());
  // // Compute directly the determinant ratio. Note: M = D^-1.
  // Scalar det_ratio = computeDetRatio(old_M, new_M);

  // EXPECT_NEAR(det_ratio, walker.getRatio(), 100 * std::numeric_limits<Scalar>::epsilon());

  // // *******************************
  // // Test vertex insertion *********
  // // *******************************
  // // Set rng values to select: first interaction vertex, tau, aux_spin,  accept
  // rng.setNewValues(std::vector<double>{0, 0.4, 0.51, 1e-6});
  // old_M = walker.getM();
  // result = walker.tryVertexInsert();
  // new_M = walker.getM();
  // ASSERT_EQ(true, result);
  // //  ASSERT_EQ(old_M.nrCols(), new_M.nrCols() - 2);
  // // Compute directly the new M.
  // walker.setMFromConfig();
  // direct_M = walker.getM();
  // for (int s = 0; s < 2; ++s)
  //   for (int j = 0; j < new_M[s].nrCols(); j++)
  //     for (int i = 0; i < new_M[s].nrRows(); i++)
  //       EXPECT_NEAR(direct_M[s](i, j), new_M[s](i, j), 20 * std::numeric_limits<Scalar>::epsilon());
  // det_ratio = computeDetRatio(old_M, new_M);
  // EXPECT_NEAR(det_ratio, walker.getRatio(), 100 * std::numeric_limits<Scalar>::epsilon());

  // // ****************************************
  // // Test last vertex removal and insertion *
  // // ****************************************
  // rng.setNewValues(std::vector<double>(100, 0));
  // const int n_vertices = walker.order();
  // for (int i = 0; i < n_vertices - 1; i++)
  //   walker.tryVertexRemoval();
  // ASSERT_EQ(1, walker.order());
  // old_M = walker.getM();
  // result = walker.tryVertexRemoval();
  // // walker.getM() is now empty
  // det_ratio = determinant(old_M) / 1.;
  // EXPECT_NEAR(det_ratio, walker.getRatio(), 1e-5);
  // // Test insertion.
  // rng.setNewValues(std::vector<double>{0, 0.5, 0, 1e-6});
  // result = walker.tryVertexInsert();
  // new_M = walker.getM();
  // det_ratio = 1. / determinant(new_M);
  // EXPECT_NEAR(det_ratio, walker.getRatio(), 1e-5);
  // walker.setMFromConfig();
  // direct_M = walker.getM();
  // for (int s = 0; s < 2; ++s)
  //   for (int j = 0; j < new_M[s].nrCols(); j++)
  //     for (int i = 0; i < new_M[s].nrRows(); i++)
  //       EXPECT_NEAR(direct_M[s](i, j), new_M[s](i, j), 10 * std::numeric_limits<Scalar>::epsilon());
}
