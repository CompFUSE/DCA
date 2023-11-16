// Copyright (C) 2023 ETH Zurich
// Copyright (C) 2023 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Wrapper to access protected members of the CT-AUX walker inside the tests.
// based on ctint/walker_wrapper.hpp

#ifndef TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_WALKER_WRAPPER_HPP
#define TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_WALKER_WRAPPER_HPP

#include "dca/distribution/dist_types.hpp"
#include "dca/linalg/device_type.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/ctaux_walker.hpp"
#include "dca/phys/dca_data/dca_data.hpp"

namespace testing {
namespace phys {
namespace solver {
namespace ctaux {
// testing::phys::solver::ctint::

using dca::DistType;

using namespace dca::phys::solver::ctaux;
template <typename SCALAR, dca::linalg::DeviceType DEVICE, class Parameters>
class CTAUXWalkerWrapper : public CtauxWalker<DEVICE, Parameters, dca::phys::DcaData<Parameters>> {
public:
  using Base = CtauxWalker<DEVICE, Parameters, dca::phys::DcaData<Parameters>>;
  using Scalar = SCALAR;
  using Real = dca::util::RealAlias<Scalar>;
  using Rng = typename Base::Rng;

  CTAUXWalkerWrapper(Parameters& parameters_ref, dca::phys::DcaData<Parameters>& data, Rng& rng_ref,
                     int id)
      : Base(parameters_ref, data, rng_ref, 0) {
    Base::initialize(0);
  }

  using Base::doStep;

  void doStep(int& some_int) {
    Base::doStep(some_int);
  }

  using Matrix = dca::linalg::Matrix<Scalar, DEVICE>;

  auto& getNUp() {
    return Base::N_up;
  }

  auto& getNDn() {
    return Base::N_dn;
  }

  auto& getGUp() {
    return Base::G_up;
  }

  auto& getGDn() {
    return Base::G_dn;
  }

  auto& getG0Up() {
    return Base::G0_up;
  }

  auto& getG0Dn() {
    return Base::G0_dn;
  }
  
  auto& getGamma_up_CPU() {
    return Base::Gamma_up_CPU;
  }

  auto& getGamma_dn_CPU() {
    return Base::Gamma_dn_CPU;
  }

  Real getRatio() const {
    return Base::det_ratio_[0] * Base::det_ratio_[1];
  }

  auto getAcceptanceProbability() const {
    return Base::acceptance_prob_;
  }

  const auto& getWalkerConfiguration() const {
    return Base::configuration_;
  }

  auto& getWarmUpdateExpansionOrder() {
    return Base::warm_up_expansion_order_;
  }

  auto& getExpVCPU() {
    return Base::exp_V_CPU;
  }

  auto& getExpV() {
    return Base::exp_V;
  }

  auto& getExpDeltaVCPU() {
    return Base::exp_delta_V_CPU;
  }

  auto& getExpDeltaV() {
    return Base::exp_delta_V;
  }

  auto& getVertexInd() {
    return Base::vertex_indixes;
  }
private:
  using Base::warm_up_expansion_order_;
  using Base::configuration_;
};

}  // namespace ctaux
}  // namespace solver
}  // namespace phys
}  // namespace testing

#endif  //  TEST_UNIT_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_WALKER_WRAPPER_HPP
