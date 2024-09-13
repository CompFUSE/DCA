// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Jérémie Bouquet   (bouquetj@gmail.com).
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch).
//
// This class tests the CPU walker used by the ctint cluster solver. The fast updated matrix
// are compared with their direct computation.

#include "dca/platform/dca_gpu.h"
using Scalar = double;
#include "test/mock_mcconfig.hpp"
namespace dca {
namespace config {
using McOptions = MockMcOptions<Scalar>;
}  // namespace config
}  // namespace dca

#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_cpu_submatrix.hpp"
#include "gtest/gtest.h"

#include "walker_wrapper.hpp"
#include "walker_wrapper_submatrix.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/details/solver_methods.hpp"

constexpr char input_name[] =
    DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/ctint/walker/submatrix_input.json";

template <typename Scalar>
using CtintWalkerSubmatrixTest =
  typename dca::testing::G0Setup<Scalar, dca::testing::LatticeBilayer, dca::ClusterSolverId::CT_INT, input_name>;

using namespace dca::phys::solver;

// Currently testing float isn't really possible due to the way the Scalar type is
// carried through from mc_options. See test_setup.hpp PD
using ScalarTypes = ::testing::Types<double>; //double, 
TYPED_TEST_CASE(CtintWalkerSubmatrixTest, ScalarTypes);

// Compare the submatrix update with a direct computation of the M matrix, and compare the
// acceptance probability to
// the CTINT walker with no submatrix update.
TYPED_TEST(CtintWalkerSubmatrixTest, doSteps) {
  using Scalar = TypeParam;
  using Parameters = typename TestFixture::Parameters;

  using Walker = testing::phys::solver::ctint::WalkerWrapper<Scalar, Parameters>;
  using Matrix = typename Walker::Matrix;
  using MatrixPair = std::array<Matrix, 2>;
  using SubmatrixWalker =
    testing::phys::solver::ctint::WalkerWrapperSubmatrix<Scalar, Parameters, dca::linalg::CPU>;

  std::vector<double> setup_rngs{0., 0.00, 0.9,  0.5, 0.01, 0,    0.75, 0.02,
                                 0,  0.6,  0.03, 1,   0.99, 0.04, 0.99};
  typename TestFixture::RngType rng(setup_rngs);

  auto& data = *TestFixture::data_;
  auto& parameters = TestFixture::parameters_;

  G0Interpolation<dca::linalg::CPU, Scalar> g0(
      dca::phys::solver::ctint::details::shrinkG0(data.G0_r_t));
  typename TestFixture::LabelDomain label_dmn;
  Walker::setDMatrixBuilder(g0);
  Walker::setDMatrixAlpha(parameters.getAlphas(), false);
  Walker::setInteractionVertices(data, parameters);

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

  for (int steps = 1; steps <= 8; ++steps) {
    rng.setNewValues(setup_rngs);
    SubmatrixWalker walker(parameters, rng);

    MatrixPair old_M(walker.getM());
    rng.setNewValues(rng_vals);
    walker.doStep(steps);
    MatrixPair new_M(walker.getM());
    auto config = walker.getWalkerConfiguration();

    // Compute directly the new M.
    walker.setMFromConfig();
    MatrixPair direct_M(walker.getM());

    using dca::util::RealAlias;
    
    // This should just be the RealAliases of the scalar since in the
    // complex case we are checking the Euclidean norm.
    const auto tolerance = 1000.0 * std::numeric_limits<RealAlias<Scalar>>::epsilon();

    for (int s = 0; s < 2; ++s)
	EXPECT_TRUE(dca::linalg::matrixop::areNear(direct_M[s], new_M[s], tolerance));

    // Compare with non submatrix walker.
    rng.setNewValues(setup_rngs);
    Walker walker_nosub(parameters, rng);

    rng.setNewValues(rng_vals);
    for (int i = 0; i < steps; ++i)
      walker_nosub.doStep();

    // this needs to be std::abs because it could be a "complex" probability
    EXPECT_NEAR(std::abs(walker.getAcceptanceProbability()), std::abs(walker_nosub.getAcceptanceProbability()),
                tolerance);

    auto config_nosubm = walker_nosub.getWalkerConfiguration();
    ASSERT_EQ(config.size(), config_nosubm.size());
    // The final configuration is the same up to a permutation.
    for (int i = 0; i < config.size(); ++i) {
      bool found = false;
      for (int j = 0; j < config_nosubm.size(); ++j)
        if (config[i] == config_nosubm[j]) {
          found = true;
          break;
        }
      EXPECT_TRUE(found);
    }
  }
}
