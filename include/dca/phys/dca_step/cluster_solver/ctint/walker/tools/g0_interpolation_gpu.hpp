// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class organizes the interpolation G0(tau) for tau in [0, beta]
// specialization for CPU.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_G0_INTERPOLATION_GPU_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_G0_INTERPOLATION_GPU_HPP
#ifdef DCA_HAVE_CUDA

#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/device_interpolation_data.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/g0_interpolation.hpp"

#include <stdexcept>

#include "dca/linalg/device_type.hpp"
#include "dca/linalg/vector.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/device_helper/ctint_helper.cuh"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/g0_interpolation.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/kernels_interface.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

template <typename Real>
class G0Interpolation<linalg::GPU, Real> final : public DeviceInterpolationData<Real>,
                                                 public G0Interpolation<linalg::CPU, Real> {
private:
  using PDmn = G0ParametersDomain;
  using PDmn0 = func::dmn_0<G0ParametersDomain>;
  using TDmn = func::dmn_0<domains::time_domain>;
  using HostInterpolation = G0Interpolation<linalg::CPU, Real>;

public:
  G0Interpolation() = default;
  G0Interpolation(const G0Interpolation& other) = delete;

  // See this->initialize(G0_pars_t).
  template <class InputDmn>
  G0Interpolation(const func::function<double, InputDmn>& G0_pars_t);
  // In: G0_pars_t. Assumed to be a function of discrete labels and time (in this order) and to
  //             be antiperiodic in time.
  template <class InputDmn>
  void initialize(const func::function<double, InputDmn>& G0_pars_t);
  // Constructor. Reshape the G0Dmn and calls the first constructor.
  // In: G0_r_t: function of dmn_variadic<D1, D2, ..., TDmn>

  // Returns cubic interpolation of G0(tau) in the spin-band-position defined by lindex.
  // Call from the CPU only for testing purposes.
  Real operator()(Real tau, int lindex) const;

  //  int getStride() const {
  //    if (time_slices_ == -1)
  //      throw(std::logic_error("The interpolation was not initialized."));
  //    return time_slices_ * COEFF_SIZE;
  //  }

  using HostInterpolation::COEFF_SIZE;

private:
  using typename HostInterpolation::CoeffDmn;
  using typename HostInterpolation::PTime0;
  using typename HostInterpolation::InterpolationDmn;

  linalg::Vector<Real, linalg::GPU> G0_coeff_;
  linalg::Vector<Real, linalg::GPU> g0_minus_dev_;
  int time_slices_ = -1;
};

template <typename Real>
template <class InputDmn>
void G0Interpolation<linalg::GPU, Real>::initialize(const func::function<double, InputDmn>& G0_pars_t) {
  HostInterpolation::initialize(G0_pars_t);
  assert(HostInterpolation::beta_);

  time_slices_ = this->getTimeSlices();
  DeviceInterpolationData<Real>::beta_ = HostInterpolation::beta_;
  DeviceInterpolationData<Real>::n_div_beta_ = HostInterpolation::n_div_beta_;
  DeviceInterpolationData<Real>::stride_ = HostInterpolation::getStride();

  G0_coeff_.resizeNoCopy(HostInterpolation::G0_coeff_.size());
  g0_minus_dev_.setAsync(HostInterpolation::g0_minus_, 0, 0);

  linalg::Vector<Real, linalg::CPU> host_coeff(HostInterpolation::G0_coeff_.size());
  for (std::size_t i = 0; i < host_coeff.size(); ++i) {
    host_coeff[i] = HostInterpolation::G0_coeff_(i);
  }
  G0_coeff_.setAsync(host_coeff, 0, 0);

  linalg::util::syncStream(0, 0);

  // Copy pointer to the data structure.
  DeviceInterpolationData<Real>::values_ = G0_coeff_.ptr();
  DeviceInterpolationData<Real>::g0_minus_ = g0_minus_dev_.ptr();
}

template <typename Real>
template <class InputDmn>
G0Interpolation<linalg::GPU, Real>::G0Interpolation(const func::function<double, InputDmn>& G0_pars_t) {
  initialize(G0_pars_t);
}

template <typename Real>
Real G0Interpolation<linalg::GPU, Real>::operator()(Real tau, int lindex) const {
  return details::interpolateSlow(tau, lindex, static_cast<DeviceInterpolationData<Real>>(*this));
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_G0_INTERPOLATION_GPU_HPP
