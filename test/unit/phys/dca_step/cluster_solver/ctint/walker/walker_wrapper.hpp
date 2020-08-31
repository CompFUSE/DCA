// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Wrapper to access protected members of the CT-INT walker inside the tests.

#ifndef TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_WALKER_WRAPPER_HPP
#define TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_WALKER_WRAPPER_HPP

#include "dca/linalg/device_type.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/ctint_walker_cpu.hpp"
#include "dca/phys/dca_data/dca_data.hpp"

namespace testing {
namespace phys {
namespace solver {
namespace ctint {
// testing::phys::solver::ctint::

using namespace dca::phys::solver::ctint;
template <class Parameters, typename Real = double>
class WalkerWrapper : public CtintWalker<dca::linalg::CPU, Parameters, Real> {
public:
  using BaseClass = CtintWalker<dca::linalg::CPU, Parameters, Real>;
  using Rng = typename BaseClass::Rng;

  WalkerWrapper(Parameters& parameters_ref, Rng& rng_ref)
      : BaseClass(parameters_ref, dca::phys::DcaData<Parameters>(parameters_ref), rng_ref, 0) {
    BaseClass::initialize(0);
  }

  using BaseClass::doStep;

  bool tryVertexInsert() {
    BaseClass::initializeStep();
    return BaseClass::tryVertexInsert();
  }
  bool tryVertexRemoval() {
    BaseClass::initializeStep();
    return BaseClass::tryVertexRemoval();
  }

  using BaseClass::setMFromConfig;
  using BaseClass::getM;

  using Matrix = dca::linalg::Matrix<Real, dca::linalg::CPU>;

  void setM(const Matrix& m) {
    BaseClass::getM() = m;
  }

  Real getRatio() const {
    return BaseClass::det_ratio_[0] * BaseClass::det_ratio_[1];
  }

  Real getAcceptanceProbability() const {
    return BaseClass::acceptance_prob_;
  }

  const auto& getWalkerConfiguration() const {
    return BaseClass::configuration_;
  }

private:
  using BaseClass::configuration_;
};

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace testing

#endif  //  TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_WALKER_WRAPPER_HPP
