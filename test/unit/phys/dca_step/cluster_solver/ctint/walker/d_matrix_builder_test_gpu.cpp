#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/d_matrix_builder_gpu.hpp"

#include "gtest/gtest.h"

#include "dca/linalg/matrixop.hpp"
#include "dca/linalg/util/cuda_stream.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/details/solver_methods.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/structs/ct_int_configuration_gpu.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

using G0Setup = dca::testing::G0Setup<dca::testing::LatticeSquare, dca::phys::solver::CT_INT>;
using namespace dca::phys::solver;
using Config = ctint::SolverConfiguration<dca::linalg::GPU>;
using HostMatrix = dca::linalg::Matrix<double, dca::linalg::CPU>;
using DeviceMatrix = dca::linalg::Matrix<double, dca::linalg::GPU>;

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

  ctint::DMatrixBuilder<dca::linalg::GPU> builder(g0, RDmn::parameter_type::get_subtract_matrix(),
                                                  label_dmn.get_branch_domain_steps(),
                                                  parameters.getAlphas());
  ctint::DMatrixBuilder<dca::linalg::CPU> builder_cpu(
      g0.get_host_interpolation(), RDmn::parameter_type::get_subtract_matrix(),
      label_dmn.get_branch_domain_steps(), parameters.getAlphas());

  HostMatrix G0;
  DeviceMatrix G0_dev;
  dca::linalg::util::CudaStream stream;

  const std::vector<int> sizes{1, 3, 8};
  const std::vector<int> deltas{1, 2, 3};
  int s(0);
  bool right_sector = true;
  for (int size : sizes) {
    right_sector = !right_sector;

    // Setup the configuration.
    Config configuration(parameters.get_beta(), BDmn::dmn_size(), G0Setup::interaction_vertices);
    for (int i = 0; i < size; i++)
      configuration.insertRandom(rng);
    configuration.upload(0);

    for (int delta : deltas) {
      s = !s;
      if (delta > size)
        continue;

      const int n_init = size - delta;
      const std::pair<int, int> matrix_size =
          right_sector ? std::make_pair(size, delta) : std::make_pair(delta, n_init);
      G0.resizeNoCopy(matrix_size);
      G0_dev.resizeNoCopy(matrix_size);

      builder_cpu.computeG0(G0, configuration.getSector(s), n_init, size, right_sector);
      builder.computeG0(G0_dev, configuration.getDeviceData(s), n_init, right_sector, stream);
      cudaStreamSynchronize(stream);

      HostMatrix G0_dev_copy(G0_dev);
      EXPECT_TRUE(dca::linalg::matrixop::areNear(G0, G0_dev_copy, 1e-10));
    }
  }
}
