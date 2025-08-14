// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch).
//
// This class tests the CPU walker used by the ctint cluster solver. The fast updated matrix
// are compared with their direct computation.

#include "dca/platform/dca_gpu.h"
#include <memory>
using Scalar = double;
#include "test/mock_mcconfig.hpp"
namespace dca {
namespace config {
using McOptions = MockMcOptions<Scalar>;
}  // namespace config
}  // namespace dca

#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_cpu_submatrix.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_cpu.hpp"
#include "dca/testing/gtest_h_w_warning_blocking.h"

#include "walker_wrapper.hpp"
#include "walker_wrapper_submatrix.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/details/solver_methods.hpp"

// This file MUST have

//   -"double-update-probability": 0
//   -"initial-configuration-size": 0
//   or walker wrapper will crash the test on construction
constexpr char input_name[] =
    DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/ctint/walker/multinsert_input.json";

template <typename SCALAR>
struct CtintWalkerMultiInsertTestT : public ::testing::Test {
  using G0Setup = dca::testing::G0SetupFromParam<SCALAR, dca::testing::LatticeBilayer,
                                            dca::ClusterSolverId::CT_INT, input_name>;
  using Parameters = typename G0Setup::Parameters;

  virtual void SetUp() {}

  void setParameters(Parameters& param_single, Parameters& param_double)
  {
    single_setup = std::make_unique<G0Setup>(param_single);
    single_setup->setUp();
    double_setup = std::make_unique<G0Setup>(param_double);
    double_setup->setUp();
  }

  virtual void TearDown() {
    single_setup.reset();
    double_setup.reset();
  }

  std::unique_ptr<G0Setup> single_setup;
  std::unique_ptr<G0Setup> double_setup;
};

template <typename Scalar>
using CtintWalkerMultiInsertTest = CtintWalkerMultiInsertTestT<Scalar>;

using Scalar = double;

using namespace dca::phys::solver;

constexpr int bands = dca::testing::LatticeHund::BANDS;

auto getVertexRng = [](int id) {
  const double n_vertices = 12.;
  return (n_vertices - id - 0.5) / n_vertices;
};

using CDA = dca::phys::ClusterDomainAliases<dca::testing::LatticeHund::DIMENSION>;
using RDmn = typename CDA::RClusterDmn;

using ScalarTypes = ::testing::Types<double>;
TYPED_TEST_CASE(CtintWalkerMultiInsertTest, ScalarTypes);

// Compare single versus double update without a submatrix update.
TYPED_TEST(CtintWalkerMultiInsertTest, NoSubmatrix) {
  using Scalar = TypeParam;
  using G0Setup = CtintWalkerMultiInsertTest<Scalar>::G0Setup;
  using Concurrency = typename G0Setup::Concurrency;
  using Parameters = typename G0Setup::Parameters;
  using Walker = typename testing::phys::solver::ctint::WalkerWrapper<Scalar, Parameters>;
  using Matrix = typename Walker::Matrix;
  using MatrixPair = std::array<Matrix, 2>;

  Concurrency concurrency(0,nullptr);

  Parameters param_single("", concurrency);
  param_single.template readInput<dca::io::JSONReader>(input_name);
  param_single.setDoubleUpdateProbability(0);
  param_single.broadcast();

  Parameters param_double("", concurrency);
  param_double.template readInput<dca::io::JSONReader>(input_name);
  param_double.setDoubleUpdateProbability(1);
  param_double.broadcast();

  this->setParameters(param_single, param_double);

  typename G0Setup::RngType rng(std::vector<Scalar>{});
  G0Interpolation<dca::linalg::CPU, Scalar> g0(
      dca::phys::solver::ctint::details::shrinkG0(this->single_setup->data_->G0_r_t));
  typename G0Setup::LabelDomain label_dmn;

  using DMatrixBuilder = dca::phys::solver::ctint::DMatrixBuilder<dca::linalg::CPU, Scalar>;
  DMatrixBuilder d_matrix_builder_single(g0, bands, RDmn());
  d_matrix_builder_single.setAlphas(param_single.getAlphas(), false);

  auto& data_single = this->single_setup->data_;

  Walker::setInteractionVertices(*data_single, param_single);

  Walker walker_single(param_single, *data_single, rng, d_matrix_builder_single);

  // interaction, tau, aux, accept
  rng.setNewValues(std::vector<double>{getVertexRng(6), 0.41, 0, 0});
  walker_single.tryVertexInsert();
  rng.setNewValues(std::vector<double>{getVertexRng(8), 0.53, 1, 0});
  walker_single.tryVertexInsert();

  // interaction, tau, aux, accept
  rng.setNewValues(std::vector<double>{getVertexRng(7), 0.34, 0, 0});
  walker_single.tryVertexInsert();
  rng.setNewValues(std::vector<double>{getVertexRng(9), 0.36, 0, 0});
  walker_single.tryVertexInsert();

  // interaction, tau, aux, accept
  rng.setNewValues(std::vector<double>{getVertexRng(6), 0.24, 1, 0});
  walker_single.tryVertexInsert();
  rng.setNewValues(std::vector<double>{getVertexRng(8), 0.19, 1, 0});
  walker_single.tryVertexInsert();

  // interaction, accept
  rng.setNewValues(std::vector<double>{0, 0});
  walker_single.tryVertexRemoval();
  rng.setNewValues(std::vector<double>{0, 0});
  walker_single.tryVertexRemoval();

  auto M1 = walker_single.getM();
  walker_single.setMFromConfig();
  auto M1_dir = walker_single.getM();
  int final_size = M1_dir[0].nrCols();
  ASSERT_EQ(M1[0].nrCols(), final_size);
  ASSERT_EQ(M1[1].nrCols(), final_size);

  const double prob_single = walker_single.getAcceptanceProbability();

  auto& data_double = this->double_setup->data_;

  Walker::setInteractionVertices(*data_double, param_double);

  DMatrixBuilder d_matrix_builder_double(g0, bands, RDmn());
  d_matrix_builder_double.setAlphas(param_double.getAlphas(), false);

  Walker walker_double(param_double, *data_double, rng, d_matrix_builder_double);
  ////////////////////////////////
  // interaction, partner, tau, aux, tau, aux, accept
  rng.setNewValues(std::vector<double>{getVertexRng(6), 0, 0.41, 0, 0.53, 1, 0});
  walker_double.tryVertexInsert();
  // interaction, partner tau, aux, tau, aux, accept
  rng.setNewValues(std::vector<double>{getVertexRng(7), 0.8, 0.34, 0, 0.36, 0, 0});
  walker_double.tryVertexInsert();
  rng.setNewValues(std::vector<double>{getVertexRng(6), 0, 0.24, 1, 0.19, 1, 0});
  walker_double.tryVertexInsert();
  // //
  // // first_id, double removal, pa rtner_id, accept
  rng.setNewValues(std::vector<double>{0, 0., 0, 0., 0.});
  walker_double.tryVertexRemoval();

  auto M2 = walker_double.getM();
  walker_double.setMFromConfig();
  auto M2_dir = walker_double.getM();
  final_size = M2_dir[0].nrCols();
  ASSERT_EQ(M2[0].nrCols(), final_size);
  ASSERT_EQ(M2[1].nrCols(), final_size);

  for (int s = 0; s < 2; ++s) {
    EXPECT_TRUE(dca::linalg::matrixop::areNear(M1[s], M1_dir[s], 1e-7));
    EXPECT_TRUE(dca::linalg::matrixop::areNear(M2[s], M2_dir[s], 1e-7));

    ASSERT_EQ(M1[s].size(), M2[s].size());
    std::vector<double> el1, el2;
    for (int j = 0; j < M1[s].nrRows(); ++j)
      for (int i = 0; i < M1[s].nrCols(); ++i) {
        el1.push_back(M1[s](i, j));
        el2.push_back(M2[s](i, j));
      }
    std::sort(el1.begin(), el1.end());
    std::sort(el2.begin(), el2.end());

    const double prob_double = walker_single.getAcceptanceProbability();
    EXPECT_NEAR(prob_single, prob_double, 1e-7);
    for (int i = 0; i < el1.size(); ++i)
      EXPECT_NEAR(el1[i], el2[i], 1e-5);
  }
}

// Compare the double insertion probability of a submatrix walker with a non-submatrix one.
// TEST_F(G0Setup, Submatrix) {
//   using Scalar = TypeParam;
//   using G0Setup = decltype(*this)::G0Setup;
//   using WalkerSubmatrix =
//     testing::phys::solver::ctint::WalkerWrapperSubmatrix<Scalar, G0Setup::Parameters>;
//   typename G0Setup::RngType rng(std::vector<double>{});
//   G0Interpolation<dca::linalg::CPU, double> g0(
//       dca::phys::solver::ctint::details::shrinkG0(this->single_setup->data_->G0_r_t);
//   typename G0Setup::LabelDomain label_dmn;

//   using DMatrixBuilder = dca::phys::solver::ctint::DMatrixBuilder<dca::linalg::CPU, Scalar>;
//   DMatrixBuilder d_matrix_builder(g0, bands, RDmn());
//   d_matrix_builder.setAlphas(parameters_.getAlphas(), false);
//   parameters_.setDoubleUpdateProbability(1);
//   Walker::setInteractionVertices(*data_, parameters_);

//   Walker walker(parameters_, rng, d_matrix_builder);

//   // interaction, partner, tau, aux, tau, aux, accept
//   rng.setNewValues(std::vector<double>{getVertexRng(6), 0, 0.41, 0, 0.53, 1, 0});
//   walker.tryVertexInsert();
//   rng.setNewValues(std::vector<double>{getVertexRng(7), 0.9, 0.34, 0, 0.36, 0, 0});
//   walker.tryVertexInsert();
//   // first idx, double rem, second partner, accept
//   rng.setNewValues(std::vector<double>{0, 0, 0, 0});
//   walker.tryVertexRemoval();

//   auto M1 = walker.getM();
//   const double prob1 = walker.getAcceptanceProbability();

//   /////////////////////////////
//   using DMatrixBuilder = dca::phys::solver::ctint::DMatrixBuilder<dca::linalg::CPU, Scalar>;
//   DMatrixBuilder d_matrix_cpu(g0, bands, RDmn());
//   d_matrix_builder.setAlphas(parameters_.getAlphas(), false);
//   WalkerSubmatrix::setInteractionVertices(*data_, parameters_);

//   WalkerSubmatrix walker_subm(parameters_, rng, d_matrix_builder);

//   std::vector<double> random_vals{// (insert, first_id, partner_id, tau, aux, tau, aux, accept) x 2
//                                   0, getVertexRng(6), 0, 0.41, 0, 0.53, 1, 0, 0, getVertexRng(7),
//                                   0.9, 0.34, 0, 0.36, 0, 0,
//                                   // with remov, first_id, double_removal, second_id, accept
//                                   1, 0, 0, 0, 0};

//   rng.setNewValues(random_vals);
//   walker_subm.doStep(3);

//   auto M2 = walker_subm.getM();
//   const double prob2 = walker_subm.getAcceptanceProbability();

//   EXPECT_NEAR(prob1, prob2, 1e-7);
//   for (int s = 0; s < 2; ++s) {
//     EXPECT_EQ(M1[s].size(), M2[s].size());
//     EXPECT_TRUE(dca::linalg::matrixop::areNear(M1[s], M2[s], 1e-7));
//     M1[s].print();
//     M2[s].print();
//   }
// }
