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

template <>
class G0Interpolation<linalg::GPU> final : public DeviceInterpolationData,
                                           public G0Interpolation<linalg::CPU> {
private:
  using PDmn = G0ParametersDomain;
  using PDmn0 = func::dmn_0<G0ParametersDomain>;
  using TDmn = func::dmn_0<domains::time_domain>;
  using HostInterpolation = G0Interpolation<linalg::CPU>;

public:
  G0Interpolation() = default;
  G0Interpolation(const G0Interpolation<linalg::GPU>& other) = delete;

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
  double operator()(double tau, int lindex) const;

  //  int getStride() const {
  //    if (time_slices_ == -1)
  //      throw(std::logic_error("The interpolation was not initialized."));
  //    return time_slices_ * COEFF_SIZE;
  //  }

  using HostInterpolation::COEFF_SIZE;

private:
  using HostInterpolation::CoeffDmn;
  using HostInterpolation::PTime0;
  using HostInterpolation::InterpolationDmn;

  linalg::Vector<double, linalg::GPU> G0_coeff_;
  linalg::Vector<double, linalg::GPU> g0_minus_dev_;
  int time_slices_ = -1;
};

template <class InputDmn>
void G0Interpolation<linalg::GPU>::initialize(const func::function<double, InputDmn>& G0_pars_t) {
  HostInterpolation::initialize(G0_pars_t);
  assert(HostInterpolation::beta_);

  time_slices_ = getTimeSlices();
  DeviceInterpolationData::beta_ = HostInterpolation::beta_;
  DeviceInterpolationData::n_div_beta_ = HostInterpolation::n_div_beta_;
  DeviceInterpolationData::stride_ = HostInterpolation::getStride();

  g0_minus_dev_.set(HostInterpolation::g0_minus_, 0, 0);
  G0_coeff_.resizeNoCopy(HostInterpolation::G0_coeff_.size());
  cudaMemcpy(G0_coeff_.ptr(), HostInterpolation::G0_coeff_.values(),
             G0_coeff_.size() * sizeof(decltype(G0_coeff_.ptr())), cudaMemcpyHostToDevice);
  // Copy pointer to the data structure.
  DeviceInterpolationData::values_ = G0_coeff_.ptr();
  DeviceInterpolationData::g0_minus_ = g0_minus_dev_.ptr();
}

template <class InputDmn>
G0Interpolation<linalg::GPU>::G0Interpolation(const func::function<double, InputDmn>& G0_pars_t) {
  initialize(G0_pars_t);
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_G0_INTERPOLATION_GPU_HPP
