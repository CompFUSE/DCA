// Copyright (C) 2024 ETH Zurich
// Copyright (C) 2024 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Jérémie Bouquet   (bouquetj@gmail.com).
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch).
//         Peter W. Doak     (doakpw@ornl.gov)
//
// This class tests the GPU walker used by the ctint cluster solver by comparing it with the CPU
// version.

#include "dca/platform/dca_gpu.h"
using Scalar = double;
#include "test/mock_mcconfig.hpp"
namespace dca {
namespace config {
using McOptions = MockMcOptions<Scalar>;
}  // namespace config
}  // namespace dca

#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_gpu_submatrix.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_cpu_submatrix.hpp"
#include "dca/testing/gtest_h_w_warning_blocking.h"

#include "walker_wrapper.hpp"
#include "walker_wrapper_submatrix.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/details/solver_methods.hpp"
#include "dca/util/to_string.hpp"

using namespace dca::addt_str_oper;

constexpr char input_name[] =
    DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/ctint/walker/submatrix_input.json";

using dca::linalg::CPU;
using dca::linalg::GPU;

template <typename SCALAR>
struct CtINTWalkerSubmatrixGPUTestT : public ::testing::Test {
  using G0Setup = dca::testing::G0SetupBare<SCALAR, dca::testing::LatticeBilayer,
                                            dca::ClusterSolverId::CT_INT, input_name>;
  virtual void SetUp() {
    host_setup.SetUp();
    gpu_setup.SetUp();
  }

  virtual void TearDown() {}
  G0Setup host_setup;
  G0Setup gpu_setup;
};

using CDA = dca::phys::ClusterDomainAliases<dca::testing::LatticeBilayer::DIMENSION>;
using RDmn = typename CDA::RClusterDmn;

using namespace dca::phys::solver;

template <typename Scalar>
using CtintWalkerSubmatrixGpuTest = CtINTWalkerSubmatrixGPUTestT<Scalar>;

template <dca::linalg::DeviceType DEVICE>
using DMatrixBuilder = dca::phys::solver::ctint::DMatrixBuilder<DEVICE, Scalar>;

// Currently testing float isn't really possible due to the way the Scalar type is
// carried through from mc_options. See test_setup.hpp PD
using ScalarTypes = ::testing::Types<double>;  // double,
TYPED_TEST_CASE(CtintWalkerSubmatrixGpuTest, ScalarTypes);

// Compare the submatrix update with a direct computation of the M matrix, and compare the
// acceptance probability to
// the CTINT walker with no submatrix update.
TYPED_TEST(CtintWalkerSubmatrixGpuTest, doSteps) {
  using Scalar = TypeParam;
  using Parameters = typename CtintWalkerSubmatrixGpuTest<Scalar>::G0Setup::Parameters;
  using Walker = testing::phys::solver::ctint::WalkerWrapper<Scalar, Parameters>;
  using Matrix = typename Walker::Matrix;
  using MatrixPair = std::array<Matrix, 2>;
  using SbmWalkerCpu =
      testing::phys::solver::ctint::WalkerWrapperSubmatrix<Scalar, Parameters, dca::linalg::CPU>;
  using SbmWalkerGpu =
      testing::phys::solver::ctint::WalkerWrapperSubmatrix<Scalar, Parameters, dca::linalg::GPU>;

  std::vector<double> setup_rngs{0., 0.00, 0.9,  0.5, 0.01, 0,    0.75, 0.02,
                                 0,  0.6,  0.03, 1,   0.99, 0.04, 0.99};
  typename CtintWalkerSubmatrixGpuTest<Scalar>::G0Setup::RngType cpu_rng(setup_rngs);
  typename CtintWalkerSubmatrixGpuTest<Scalar>::G0Setup::RngType gpu_rng(setup_rngs);

  auto& cpu_data = *this->host_setup.data_;
  auto& gpu_data = *this->gpu_setup.data_;
  auto& cpu_parameters = this->host_setup.parameters_;
  auto& gpu_parameters = this->gpu_setup.parameters_;

  // beta of 1 can result in expansion factor of 0 with default bilayer model params
  gpu_parameters.set_beta(4);
  cpu_parameters.set_beta(4);

  const auto g0_func_cpu = dca::phys::solver::ctint::details::shrinkG0(cpu_data.G0_r_t);
  G0Interpolation<CPU, Scalar> g0_cpu(g0_func_cpu);
  const auto g0_func_gpu = dca::phys::solver::ctint::details::shrinkG0(gpu_data.G0_r_t);
  G0Interpolation<GPU, Scalar> g0_gpu(g0_func_gpu);
  typename CtintWalkerSubmatrixGpuTest<Scalar>::G0Setup::LabelDomain label_dmn;

  constexpr int bands = dca::testing::LatticeBilayer::BANDS;

  DMatrixBuilder<dca::linalg::CPU> d_matrix_cpu(g0_cpu, bands, RDmn());
  d_matrix_cpu.setAlphas(cpu_parameters.getAlphas(), false);  // cpu_parameters.adjustAlphaDd());
  SbmWalkerCpu::setInteractionVertices(cpu_data, cpu_parameters);
  DMatrixBuilder<dca::linalg::GPU> d_matrix_gpu(g0_gpu, bands, RDmn());
  d_matrix_gpu.setAlphas(gpu_parameters.getAlphas(), false);  // gpu_parameters.adjustAlphaDd());
  SbmWalkerGpu::setInteractionVertices(cpu_data, gpu_parameters);

  // ************************************
  // Test vertex insertion / removal ****
  // ************************************
  // Set rng values.
  //
  // Insertion, vertex_id, tau, aux_spin, acceptance_rng
  // Removal, vertex_id, acceptance_rng
  // ...
  // Note: if acceptance_rng <= 0 the move is always accepted, if it is > 1 the move is always
  // rejected.
  const std::vector<double> rng_vals{
      0, 0,    0.1, 0.8, -1,  // Insertion.
      0, 0.99, 0.2, 0.8, -1,  // Insertion.
      0, 0,    0.3, 0.8, 2,   // Insertion. Rejected.
      1, 0,    -1,            // Remove pre-existing.
      1, 0.99, -1,            // Remove recently inserted.
      1, 0.99, 2,             // Remove recently inserted. Rejected
      1, 0,    2,             // Remove . Rejected
      0, 0.99, 0.4, 0.2, -1,  // Insertion
  };

  for (const int initial_size : std::array<int, 1>{5}) {
    cpu_parameters.setInitialConfigurationSize(initial_size);
    gpu_parameters.setInitialConfigurationSize(initial_size);
    for (int steps = 1; steps <= 8; ++steps) {
      cpu_rng.setNewValues(setup_rngs);
      SbmWalkerCpu walker_cpu(cpu_parameters, cpu_rng, d_matrix_cpu);
      gpu_rng.setNewValues(setup_rngs);
      SbmWalkerGpu walker_gpu(gpu_parameters, gpu_rng, d_matrix_gpu);

      // I don't think we can call these before steps are done.  At least the CPU implementation has nan's in M
      // MatrixPair old_M_cpu(walker_cpu.getM());
      // MatrixPair old_M_gpu(walker_gpu.getM());

      constexpr Scalar tolerance = std::numeric_limits<Scalar>::epsilon() * 100;

      // for (int s = 0; s < 2; ++s)
      //   EXPECT_TRUE(dca::linalg::matrixop::areNear(old_M_cpu[s], old_M_gpu[s], tolerance));

      cpu_rng.setNewValues(rng_vals);
      walker_cpu.doStep(steps);
      MatrixPair new_M_cpu = walker_cpu.getM();

      gpu_rng.setNewValues(rng_vals);
      walker_gpu.uploadConfiguration();
      walker_gpu.doStep(steps);

      auto M_cpu = walker_cpu.getM();
      auto M_gpu = walker_gpu.getM();
      // for (int s = 0; s < 2; ++s)
      //   EXPECT_TRUE(dca::linalg::matrixop::areNear(M_cpu[s], M_gpu[s], tolerance));

      // The final configuration is the same.
      const auto& config1 = walker_cpu.getWalkerConfiguration();
      const auto& config2 = walker_gpu.getWalkerConfiguration();
      ASSERT_EQ(config1.size(), config2.size());
      for (int i = 0; i < config1.size(); ++i)
        EXPECT_EQ(config1[i], config2[i]);
      EXPECT_EQ(walker_cpu.get_sign(), walker_gpu.get_sign());
    }
  }
}

TYPED_TEST(CtintWalkerSubmatrixGpuTest, computeMInit) {
  using Scalar = TypeParam;
  using Parameters = typename CtintWalkerSubmatrixGpuTest<Scalar>::G0Setup::Parameters;
  using Walker = testing::phys::solver::ctint::WalkerWrapper<Scalar, Parameters>;
  using Matrix = typename Walker::Matrix;
  using MatrixPair = std::array<Matrix, 2>;
  using SbmWalkerCpu =
      testing::phys::solver::ctint::WalkerWrapperSubmatrix<Scalar, Parameters, dca::linalg::CPU>;
  using SbmWalkerGpu =
      testing::phys::solver::ctint::WalkerWrapperSubmatrix<Scalar, Parameters, dca::linalg::GPU>;

  std::vector<double> setup_rngs{0., 0.00, 0.9,  0.5, 0.01, 0,    0.75, 0.02,
                                 0,  0.6,  0.03, 1,   0.99, 0.04, 0.99};
  typename CtintWalkerSubmatrixGpuTest<Scalar>::G0Setup::RngType cpu_rng(setup_rngs);
  typename CtintWalkerSubmatrixGpuTest<Scalar>::G0Setup::RngType gpu_rng(setup_rngs);

  auto& cpu_data = *this->host_setup.data_;
  auto& gpu_data = *this->gpu_setup.data_;
  auto& cpu_parameters = this->host_setup.parameters_;
  auto& gpu_parameters = this->gpu_setup.parameters_;

  // beta of 1 can result in expansion factor of 0 with default bilayer model params
  // gpu_parameters.set_beta(4);
  // cpu_parameters.set_beta(4);
  
  const auto g0_func_cpu = dca::phys::solver::ctint::details::shrinkG0(cpu_data.G0_r_t);
  G0Interpolation<CPU, Scalar> g0_cpu(g0_func_cpu);
  const auto g0_func_gpu = dca::phys::solver::ctint::details::shrinkG0(gpu_data.G0_r_t);
  G0Interpolation<GPU, Scalar> g0_gpu(g0_func_gpu);
  typename CtintWalkerSubmatrixGpuTest<Scalar>::G0Setup::LabelDomain label_dmn;

  constexpr int bands = dca::testing::LatticeBilayer::BANDS;

  DMatrixBuilder<dca::linalg::CPU> d_matrix_cpu(g0_cpu, bands, RDmn());
  d_matrix_cpu.setAlphas(cpu_parameters.getAlphas(), false);  // cpu_parameters.adjustAlphaDd());
  SbmWalkerCpu::setInteractionVertices(cpu_data, cpu_parameters);
  DMatrixBuilder<dca::linalg::GPU> d_matrix_gpu(g0_gpu, bands, RDmn());
  d_matrix_gpu.setAlphas(gpu_parameters.getAlphas(), false);  // gpu_parameters.adjustAlphaDd());
  SbmWalkerGpu::setInteractionVertices(cpu_data, gpu_parameters);

  // ************************************
  // Test vertex insertion / removal ****
  // ************************************
  // Set rng values.
  //
  // Insertion, vertex_id, tau, aux_spin, acceptance_rng
  // Removal, vertex_id, acceptance_rng
  // ...
  // Note: if acceptance_rng <= 0 the move is always accepted, if it is > 1 the move is always
  // rejected.
  const std::vector<double> rng_vals{
      0, 0,    0.1, 0.8, -1,  // Insertion.
      0, 0.99, 0.2, 0.8, -1,  // Insertion.
      0, 0,    0.3, 0.8, 2,   // Insertion. Rejected.
      1, 0,    -1,            // Remove pre-existing.
      1, 0.99, -1,            // Remove recently inserted.
      1, 0.99, 2,             // Remove recently inserted. Rejected
      1, 0,    2,             // Remove . Rejected
      0, 0.99, 0.4, 0.2, -1,  // Insertion
  };

  int initial_size = 2;

  cpu_parameters.setInitialConfigurationSize(initial_size);
  gpu_parameters.setInitialConfigurationSize(initial_size);
  int steps = 2;

  cpu_rng.setNewValues(setup_rngs);
  SbmWalkerCpu walker_cpu(cpu_parameters, cpu_rng, d_matrix_cpu);
  gpu_rng.setNewValues(setup_rngs);
  SbmWalkerGpu walker_gpu(gpu_parameters, gpu_rng, d_matrix_gpu);

  constexpr Scalar tolerance = std::numeric_limits<Scalar>::epsilon() * 100;

  cpu_rng.setNewValues(rng_vals);
  walker_cpu.generateDelayedMoves(steps);
  walker_cpu.computeMInit();

  walker_gpu.uploadConfiguration();
  gpu_rng.setNewValues(rng_vals);
  walker_gpu.generateDelayedMoves(steps);
  walker_gpu.uploadConfiguration();
  walker_gpu.computeMInit();

  auto& dmatrix_gpu = walker_gpu.getDMatrixBuilder();
  EXPECT_TRUE(dca::phys::solver::details::SolverHelper::initialized());
  int nc;
  cudaMemcpyFromSymbol(&nc, &dca::phys::solver::details::cluster_real_helper.nc_, sizeof(int), 0);
  std::cout << "nc_ : " << nc << '\n';
  std::vector<int> add_matrix(nc * nc);
  int ** int_ptr;
  auto* ptr_crh = &dca::phys::solver::details::cluster_real_helper;
  cudaMemcpy(&int_ptr, dca::phys::solver::details::cluster_real_helper.add_matrix_, sizeof(int*), cudaMemcpyDeviceToHost);
  cudaMemcpy(add_matrix.data(), int_ptr, sizeof(int) * nc * nc, cudaMemcpyDeviceToHost);
  std::cout << "Add matrix:\n" << dca::vectorToString(add_matrix) << '\n';
  std::vector<int> sub_matrix(nc * nc);
  cudaMemcpy(&int_ptr, dca::phys::solver::details::cluster_real_helper.sub_matrix_, sizeof(int*), cudaMemcpyDeviceToHost);
  cudaMemcpy(sub_matrix.data(), int_ptr, sizeof(int) * nc * nc, cudaMemcpyDeviceToHost);
  std::cout << "Sub matrix:\n" << dca::vectorToString(sub_matrix) << '\n';

  std::cout << "SiteDiff: \n";
  dmatrix_gpu.getSiteDiff().print();
  
  //walker_gpu.computeGInit();
  auto M_cpu = walker_cpu.getRawM();
  auto M_gpu = walker_gpu.getRawM();

  std::cout << "Minit init cpu size " << M_cpu[0].size().first << "," << M_cpu[0].size().second << " : " <<
    M_cpu[1].size().first << "," << M_cpu[1].size().second << '\n';
  std::cout << "Minit init gpu size " << M_gpu[0].size().first << "," << M_gpu[0].size().second << " : " <<
    M_gpu[1].size().first << "," << M_gpu[1].size().second << '\n';

  std::cout << walker_cpu.getF()[-1][0] << "," << walker_cpu.getF()[-1][1] << '\n';
  std::cout << walker_cpu.getF()[1][0] << "," << walker_cpu.getF()[1][1] << '\n';
  std::cout << walker_gpu.getF()[-1][0] << "," << walker_gpu.getF()[-1][1] << '\n';
  std::cout << walker_gpu.getF()[1][0] << "," << walker_gpu.getF()[1][1] << '\n';
  std::cout << walker_gpu.getFValues()[0][0] << "," << walker_gpu.getFValues()[0][1] << '\n';
  std::cout << walker_gpu.getFValues()[1][0] << "," << walker_gpu.getFValues()[1][1] << '\n';

  M_cpu[0].print();
  M_cpu[1].print();
  M_gpu[0].print();
  M_gpu[1].print();

  for (int s = 0; s < 2; ++s) {
    int fail = 0;
    std::cout << "M_cpu[" << s << "] size: " << M_cpu[s].size() << '\n';
    std::cout << "M_gpu[" << s << "] size: " << M_gpu[s].size() << '\n';
    for (int j = 0; j < M_cpu[0].size().second; ++j)
      for (int i = 0; i < M_cpu[0].size().first; ++i) {
        if (std::abs(M_cpu[s](i, j) - M_gpu[s](i, j)) > tolerance) {
          std::cout << "M element: " << s << "," << i << "," << j << "  cpu(" << M_cpu[s](i, j)
                    << ") - gpu(" << M_gpu[s](i, j) << ") > " << tolerance << '\n';
          ++fail;
        }
      }
  }
  // EXPECT_TRUE(dca::linalg::matrixop::areNear(M_cpu[s], M_gpu[s], tolerance));

  // The final configuration is the same.
  const auto& config1 = walker_cpu.getWalkerConfiguration();
  const auto& config2 = walker_gpu.getWalkerConfiguration();
  ASSERT_EQ(config1.size(), config2.size());
  for (int i = 0; i < config1.size(); ++i)
    EXPECT_EQ(config1[i], config2[i]);
  EXPECT_EQ(walker_cpu.get_sign(), walker_gpu.get_sign());
}
