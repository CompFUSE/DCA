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

  constexpr Scalar tolerance = 10e-12;
  auto compareSubMatrix = [tolerance](auto& mat_cpu, auto& mat_gpu, const std::string& mat_tag) {
    int fail = 0;
    for (int s = 0; s < 2; ++s) {
      std::cout << mat_tag << "_cpu[" << s << "] size: " << mat_cpu[s].size() << '\n';
      std::cout << mat_tag << "_gpu[" << s << "] size: " << mat_gpu[s].size() << '\n';
      for (int j = 0; j < mat_cpu[s].size().second; ++j)
        for (int i = 0; i < mat_cpu[s].size().first; ++i) {
          if (std::abs(mat_cpu[s](i, j) - mat_gpu[s](i, j)) > tolerance) {
            std::cout << mat_tag << " element: " << s << "," << i << "," << j << "  cpu("
                      << mat_cpu[s](i, j) << ") - gpu(" << mat_gpu[s](i, j) << ") > " << tolerance
                      << '\n';
            ++fail;
          }
        }
    }
    return fail;
  };

  for (const int initial_size : std::array<int, 1>{2}) {
    cpu_parameters.setInitialConfigurationSize(initial_size);
    gpu_parameters.setInitialConfigurationSize(initial_size);
    for (int steps = 2; steps <= 8; ++steps) {
      cpu_rng.setNewValues(setup_rngs);
      SbmWalkerCpu walker_cpu(cpu_parameters, cpu_rng, d_matrix_cpu);
      gpu_rng.setNewValues(setup_rngs);
      SbmWalkerGpu walker_gpu(gpu_parameters, gpu_rng, d_matrix_gpu);

      cpu_rng.setNewValues(rng_vals);
      walker_cpu.doStep(steps);

      gpu_rng.setNewValues(rng_vals);
      walker_gpu.doStep(steps);

      // doSweep does this
      walker_gpu.uploadConfiguration();

      // The final configuration is the same.
      const auto& config1 = walker_cpu.getWalkerConfiguration();
      const auto& config2 = walker_gpu.getWalkerConfiguration();
      ASSERT_EQ(config1.size(), config2.size());
      for (int i = 0; i < config1.size(); ++i)
        EXPECT_EQ(config1[i], config2[i]);
      EXPECT_EQ(walker_cpu.get_sign(), walker_gpu.get_sign());

      auto M_cpu = walker_cpu.getM();
      auto M_gpu = walker_gpu.getM();

      auto fail = compareSubMatrix(M_cpu, M_gpu, "M");
      EXPECT_FALSE(fail);
    }
  }
}

TYPED_TEST(CtintWalkerSubmatrixGpuTest, fineGrained) {
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
  int steps = 4;

  cpu_rng.setNewValues(setup_rngs);
  SbmWalkerCpu walker_cpu(cpu_parameters, cpu_rng, d_matrix_cpu);
  gpu_rng.setNewValues(setup_rngs);
  SbmWalkerGpu walker_gpu(gpu_parameters, gpu_rng, d_matrix_gpu);

  constexpr Scalar tolerance = 10e-12;

  cpu_rng.setNewValues(rng_vals);
  walker_cpu.generateDelayedMoves(steps);
  walker_cpu.computeMInit();

  walker_gpu.uploadConfiguration();
  gpu_rng.setNewValues(rng_vals);
  walker_gpu.generateDelayedMoves(steps);
  walker_gpu.uploadConfiguration();
  walker_gpu.computeMInit();

  // walker_gpu.computeGInit();
  auto M_cpu = walker_cpu.getRawM();
  auto M_gpu = walker_gpu.getRawM();

  std::cout << "Minit init cpu size " << M_cpu[0].size().first << "," << M_cpu[0].size().second
            << " : " << M_cpu[1].size().first << "," << M_cpu[1].size().second << '\n';
  std::cout << "Minit init gpu size " << M_gpu[0].size().first << "," << M_gpu[0].size().second
            << " : " << M_gpu[1].size().first << "," << M_gpu[1].size().second << '\n';

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

  auto compareSubMatrix = [tolerance](auto& mat_cpu, auto& mat_gpu, const std::string& mat_tag) {
    int fail = 0;
    for (int s = 0; s < 2; ++s) {
      std::cout << mat_tag << "_cpu[" << s << "] size: " << mat_cpu[s].size() << '\n';
      std::cout << mat_tag << "_gpu[" << s << "] size: " << mat_gpu[s].size() << '\n';
      for (int j = 0; j < mat_cpu[s].size().second; ++j)
        for (int i = 0; i < mat_cpu[s].size().first; ++i) {
          if (std::abs(mat_cpu[s](i, j) - mat_gpu[s](i, j)) > tolerance) {
            std::cout << mat_tag << " element: " << s << "," << i << "," << j << "  cpu("
                      << mat_cpu[s](i, j) << ") - gpu(" << mat_gpu[s](i, j) << ") > " << tolerance
                      << '\n';
            ++fail;
          }
        }
    }
    return fail;
  };

  auto fail = compareSubMatrix(M_cpu, M_gpu, "M");
  EXPECT_FALSE(fail);

  walker_cpu.computeGInit();
  walker_gpu.computeGInit();
  walker_gpu.synchronize();

  auto G_cpu = walker_cpu.getRawG();
  auto G_gpu = walker_gpu.getRawG();

  fail = compareSubMatrix(G_cpu, G_gpu, "G");
  EXPECT_FALSE(fail);

  walker_cpu.mainSubmatrixProcess();
  walker_gpu.mainSubmatrixProcess();

  walker_cpu.updateM();
  walker_gpu.updateM();

  M_cpu = walker_cpu.getRawM();
  M_gpu = walker_gpu.getRawM();

  std::cout << "M post update cpu size " << M_cpu[0].size().first << "," << M_cpu[0].size().second
            << " : " << M_cpu[1].size().first << "," << M_cpu[1].size().second << '\n';
  std::cout << "M post update gpu size " << M_gpu[0].size().first << "," << M_gpu[0].size().second
            << " : " << M_gpu[1].size().first << "," << M_gpu[1].size().second << '\n';

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

  fail = compareSubMatrix(M_cpu, M_gpu, "updated_M");
  EXPECT_FALSE(fail);

  // The final configuration is the same.
  const auto& config1 = walker_cpu.getWalkerConfiguration();
  const auto& config2 = walker_gpu.getWalkerConfiguration();
  ASSERT_EQ(config1.size(), config2.size());
  for (int i = 0; i < config1.size(); ++i)
    EXPECT_EQ(config1[i], config2[i]);
  EXPECT_EQ(walker_cpu.get_sign(), walker_gpu.get_sign());

  walker_gpu.uploadConfiguration();

  auto comp_M_cpu = walker_cpu.getM();
  auto comp_M_gpu = walker_gpu.getM();

  std::cout << "computed M cpu size " << comp_M_cpu[0].size().first << ","
            << comp_M_cpu[0].size().second << " : " << comp_M_cpu[1].size().first << ","
            << comp_M_cpu[1].size().second << '\n';
  std::cout << "computed M gpu size " << comp_M_gpu[0].size().first << ","
            << comp_M_gpu[0].size().second << " : " << comp_M_gpu[1].size().first << ","
            << comp_M_gpu[1].size().second << '\n';

  std::cout << walker_cpu.getF()[-1][0] << "," << walker_cpu.getF()[-1][1] << '\n';
  std::cout << walker_cpu.getF()[1][0] << "," << walker_cpu.getF()[1][1] << '\n';
  std::cout << walker_gpu.getF()[-1][0] << "," << walker_gpu.getF()[-1][1] << '\n';
  std::cout << walker_gpu.getF()[1][0] << "," << walker_gpu.getF()[1][1] << '\n';
  std::cout << walker_gpu.getFValues()[0][0] << "," << walker_gpu.getFValues()[0][1] << '\n';
  std::cout << walker_gpu.getFValues()[1][0] << "," << walker_gpu.getFValues()[1][1] << '\n';

  comp_M_cpu[0].print();
  comp_M_cpu[1].print();
  comp_M_gpu[0].print();
  comp_M_gpu[1].print();

  fail = compareSubMatrix(comp_M_cpu, comp_M_gpu, "comp_M");
  EXPECT_FALSE(fail);
}
