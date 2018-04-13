#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/d_matrix_builder.hpp"
#include "gtest/gtest.h"

#include "dca/phys/dca_step/cluster_solver/ctint/details/solver_methods.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_configuration_gpu.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

using G0Setup = dca::testing::G0Setup<dca::testing::LatticeSquare, dca::phys::solver::CT_INT>;
using namespace dca::phys::solver;
using Config = ctint::SolverConfiguration<dca::linalg::GPU>;
using HostMatrixPair = std::array<dca::linalg::Matrix<double, dca::linalg::CPU>, 2>;
using DeviceMatrixPair = std::array<dca::linalg::Matrix<double, dca::linalg::GPU>, 2>;

TEST_F(G0Setup, RemoveAndInstertVertex) {
  // Setup rng values
  std::vector<double> rng_values(100);
  for (double& x : rng_values)
    x = double(std::rand()) / RAND_MAX;
  G0Setup::RngType rng(rng_values);

  using namespace dca::testing;
  dca::func::dmn_variadic<BDmn, BDmn, RDmn> label_dmn;

  // Setup interpolation and matrix builder class.
  ctint::G0Interpolation<dca::linalg::GPU> g0(
      dca::phys::solver::ctint::details::shrinkG0(data->G0_r_t));
  ctint::G0Interpolation<dca::linalg::CPU> g0_cpu(
      dca::phys::solver::ctint::details::shrinkG0(data->G0_r_t));

  ctint::DMatrixBuilder<dca::linalg::GPU> builder(g0, RDmn::parameter_type::get_subtract_matrix(),
                                                  label_dmn.get_branch_domain_steps(),
                                                  parameters.getAlphas());
  ctint::DMatrixBuilder<dca::linalg::CPU> builder_cpu(
      g0_cpu, RDmn::parameter_type::get_subtract_matrix(), label_dmn.get_branch_domain_steps(),
      parameters.getAlphas());

  HostMatrixPair Q, R, S;
  DeviceMatrixPair Q_dev, R_dev, S_dev;
  HostMatrixPair Q2, R2, S2;

  const std::vector<int> sizes{1, 3, 8};
  for (int size : sizes) {
    // Setup the configuration.
    Config configuration(parameters.get_beta(), BDmn::dmn_size(), G0Setup::interaction_vertices);

    for (int i = 0; i < size; i++)
      configuration.insertRandom(rng);

    builder.buildSQR(S_dev, Q_dev, R_dev, configuration, 0);
    builder_cpu.buildSQR(S2, Q2, R2, configuration);

    for (int s = 0; s < 2; ++s) {
      Q[s] = Q_dev[s];
      R[s] = R_dev[s];
      S[s] = S_dev[s];
      EXPECT_TRUE(dca::linalg::matrixop::areNear(S[s], S2[s], 1e-12));
      EXPECT_TRUE(dca::linalg::matrixop::areNear(Q[s], Q2[s], 1e-12));
      EXPECT_TRUE(dca::linalg::matrixop::areNear(R[s], R2[s], 1e-12));
    }
  }
}
