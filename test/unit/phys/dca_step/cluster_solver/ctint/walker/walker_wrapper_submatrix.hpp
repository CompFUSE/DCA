// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Jérémie Bouquet (bouquetj@gmail.com)
//
//

#ifndef TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_WALKER_WRAPPER_SUBMATRIX_HPP
#define TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_WALKER_WRAPPER_SUBMATRIX_HPP

#include "dca/linalg/device_type.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_cpu_submatrix.hpp"
#ifdef DCA_HAVE_CUDA
#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_gpu_submatrix.hpp"
#endif

namespace testing {
namespace phys {
namespace solver {
namespace ctint {
// testing::phys::solver::ctint::

using namespace dca::phys::solver::ctint;
using dca::linalg::CPU;
using dca::linalg::GPU;
using dca::linalg::DeviceType;

template <class Parameters, dca::linalg::DeviceType device, typename Scalar>
struct WalkerSelector;

template <class Parameters, typename Scalar>
struct WalkerSelector<Parameters, CPU, Scalar> {
  // Fix rng order for testing.
  using type = CtintWalkerSubmatrixCpu<Parameters, Scalar>;
};

#ifdef DCA_HAVE_CUDA
template <class Parameters, typename Scalar>
struct WalkerSelector<Parameters, GPU, Scalar> {
  using type = CtintWalkerSubmatrixGpu<Parameters, Scalar>;
};
#endif  // DCA_HAVE_CUDA

using namespace dca::phys::solver::ctint;
template <class Parameters, typename Scalar = double, DeviceType device_t = CPU>
struct WalkerWrapperSubmatrix : public WalkerSelector<Parameters, device_t, Scalar>::type {
  using BaseClass = typename WalkerSelector<Parameters, device_t, Scalar>::type;
  using typename BaseClass::Rng;
  using typename BaseClass::Data;

  WalkerWrapperSubmatrix(/*const*/ Parameters& parameters_ref, Rng& rng_ref)
      : BaseClass(parameters_ref, dca::phys::DcaData<Parameters>(parameters_ref), rng_ref, 0),
        streams_(3) {
    BaseClass::initialize(0);
  }

  void doStep(const int n_steps_to_delay) {
    BaseClass::doStep(n_steps_to_delay);
  }

  using Matrix = dca::linalg::Matrix<Scalar, CPU>;
  using MatrixPair = std::array<Matrix, 2>;

  MatrixPair getM() {
    std::array<dca::linalg::Matrix<Scalar, device_t>, 2> M;

    BaseClass::computeM(M);
#ifdef DCA_HAVE_CUDA
    cudaDeviceSynchronize();
#endif

    std::array<dca::linalg::Matrix<Scalar, CPU>, 2> M_copy{M[0], M[1]};
    return M_copy;
  }

  using BaseClass::setMFromConfig;

  const auto& getWalkerConfiguration() const {
    return BaseClass::configuration_;
  }

  auto getAcceptanceProbability() const {
    return BaseClass::acceptance_prob_;
  }

private:
  std::vector<dca::linalg::util::CudaStream> streams_;
};

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace testing

#endif  //  TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_WALKER_WRAPPER_SUBMATRIX_HPP
