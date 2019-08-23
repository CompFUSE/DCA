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
template <class Parameters, dca::linalg::DeviceType device_t = dca::linalg::CPU>
struct WalkerWrapperSubmatrix : public CtintWalkerSubmatrix<device_t, Parameters> {
  using BaseClass = CtintWalkerSubmatrix<device_t, Parameters>;
  using Rng = typename BaseClass::Rng;
  using Data = typename BaseClass::Data;

  WalkerWrapperSubmatrix(/*const*/ Parameters& parameters_ref, Rng& rng_ref)
      : BaseClass(parameters_ref, dca::phys::DcaData<Parameters>(parameters_ref), rng_ref, 0) {
      //cudaDeviceSynchronize();
      BaseClass::initialize();
      //cudaDeviceSynchronize();
  }

  void doStep(const int n_steps_to_delay) {
    BaseClass::doStep(n_steps_to_delay);
  }

  using Matrix = dca::linalg::Matrix<double, dca::linalg::CPU>;
  using MatrixPair = std::array<Matrix, 2>;

  MatrixPair getM() {
    std::array<dca::linalg::Matrix<double, device_t>, 2> M;
    std::vector<dca::linalg::util::CudaStream> s(2);
    std::vector<dca::linalg::util::CudaStream*> s_ptr{&s[0], &s[1]};
    cudaDeviceSynchronize();
    BaseClass::computeM(M, s_ptr);
    cudaDeviceSynchronize();

    std::array<dca::linalg::Matrix<double, dca::linalg::CPU>, 2> M_copy{M[0], M[1]};
    return M_copy;
  }

  void setMFromConfig() {
    BaseClass::setMFromConfig();
    BaseClass::transformM();
  }

  const auto& getWalkerConfiguration() const {
    return BaseClass::configuration_;
  }

  double getAcceptanceProbability() const {
    return BaseClass::acceptance_prob_;
  }
};

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace testing

#endif  //  TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_WALKER_WRAPPER_SUBMATRIX_HPP
