#include <dca/phys/dca_step/cluster_solver/ctint/structs/device_configuration_manager.hpp>
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/d_matrix_builder_gpu.hpp"

#include "gtest/gtest.h"

#include "dca/linalg/matrixop.hpp"
#include "dca/linalg/util/cuda_stream.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/details/solver_methods.hpp"
#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

template <typename Real>
using DMatrixBuilderGpuTest =
    dca::testing::G0Setup<dca::testing::LatticeSquare, dca::phys::solver::CT_INT>;
using namespace dca::phys::solver;

using dca::linalg::Matrix;
using dca::linalg::CPU;
using dca::linalg::GPU;

using FloatingPointTypes = ::testing::Types<float, double>;
TYPED_TEST_CASE(DMatrixBuilderGpuTest, FloatingPointTypes);

TYPED_TEST(DMatrixBuilderGpuTest, RemoveAndInstertVertex) {
  using Rng = typename TestFixture::RngType;
  using BDmn = typename TestFixture::BDmn;
  using RDmn = typename TestFixture::RDmn;
  using Real = TypeParam;

  // Setup rng values
  std::vector<double> rng_values(100);
  for (double& x : rng_values)
    x = double(std::rand()) / RAND_MAX;
  Rng rng(rng_values);

  using namespace dca::testing;
  dca::func::dmn_variadic<BDmn, BDmn, RDmn> label_dmn;

  const auto& parameters = TestFixture::parameters_;
  const auto& data = *TestFixture::data_;

  // Setup interpolation and matrix builder class.
  ctint::G0Interpolation<dca::linalg::GPU, Real> g0(
      dca::phys::solver::ctint::details::shrinkG0(data.G0_r_t));

  const int nb = BDmn::dmn_size();

  ctint::DMatrixBuilder<dca::linalg::GPU, Real> builder(g0, nb, RDmn());
  builder.setAlphas(parameters.getAlphas(), false);

  ctint::DMatrixBuilder<dca::linalg::CPU, Real> builder_cpu(g0, nb, RDmn());
  builder_cpu.setAlphas(parameters.getAlphas(), false);

  ctint::InteractionVertices interaction_vertices;
  interaction_vertices.initializeFromHamiltonian(data.H_interactions);

  Matrix<Real, CPU> G0;
  Matrix<Real, GPU> G0_dev;
  dca::linalg::util::CudaStream stream;

  const std::vector<int> sizes{1, 3, 8};
  const std::vector<int> deltas{1, 2, 3};
  int s(0);
  bool right_sector = true;
  for (int size : sizes) {
    right_sector = !right_sector;

    // Setup the configuration.
    ctint::SolverConfiguration configuration(parameters.get_beta(), BDmn::dmn_size(),
                                             interaction_vertices);
    ctint::DeviceConfigurationManager device_config;
    for (int i = 0; i < size; i++)
      configuration.insertRandom(rng);
    device_config.upload(configuration, 0);

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
      builder.computeG0(G0_dev, device_config.getDeviceData(s), n_init, right_sector, stream);
      cudaStreamSynchronize(stream);

      Matrix<Real, CPU> G0_dev_copy(G0_dev);
      constexpr Real tolerance = 100 * std::numeric_limits<Real>::epsilon();
      EXPECT_TRUE(dca::linalg::matrixop::areNear(G0, G0_dev_copy, tolerance));
    }
  }
}
