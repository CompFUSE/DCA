// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class organizes the interpolation G0(tau) for tau in [0, beta]
// specialization for CPU.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_G0_INTERPOLATION_GPU_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_G0_INTERPOLATION_GPU_HPP
#ifdef DCA_HAVE_CUDA

#include <stdexcept>

#include "dca/linalg/device_type.hpp"
#include "dca/linalg/vector.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/device_memory/global_memory_manager.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/device_interpolation_data.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/g0_interpolation.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/kernels_interface.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

template <>
class G0Interpolation<linalg::GPU> : public details::DeviceInterpolationData {
private:
  using PDmn = G0ParametersDomain;
  using PDmn0 = func::dmn_0<G0ParametersDomain>;
  using TDmn = func::dmn_0<domains::time_domain>;
  using BaseClass = details::DeviceInterpolationData;

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
  double inline operator()(double tau, int lindex) const;

  int getStride() const {
    if (time_slices_ == -1)
      throw(std::logic_error("The interpolation was not initialized."));
    return time_slices_ * COEFF_SIZE;
  }

  const auto& get_host_interpolation() const {
    return g0_cpu_;
  }

  constexpr static int COEFF_SIZE = G0Interpolation<linalg::CPU>::COEFF_SIZE;

private:
  using CoeffDmn = func::dmn_0<func::dmn<COEFF_SIZE>>;
  using PTime0 = func::dmn_0<PositiveTimeDomain>;
  using InterpolationDmn = func::dmn_variadic<CoeffDmn, PTime0, PDmn0>;

  linalg::Vector<double, linalg::GPU> G0_coeff_;
  linalg::Vector<double, linalg::GPU> g0_minus_;
  int time_slices_ = -1;

  G0Interpolation<linalg::CPU> g0_cpu_;
};

template <class InputDmn>
void G0Interpolation<linalg::GPU>::initialize(const func::function<double, InputDmn>& G0_pars_t) {
  g0_cpu_.initialize(G0_pars_t);
  assert(g0_cpu_.beta_);
  time_slices_ = g0_cpu_.getTimeSlices();
  BaseClass::beta_ = g0_cpu_.beta_;
  BaseClass::n_div_beta_ = g0_cpu_.n_div_beta_;

  g0_minus_.set(g0_cpu_.g0_zero_minus, 0, 0);
  G0_coeff_.resizeNoCopy(g0_cpu_.G0_coeff_.size());
  cudaMemcpy(G0_coeff_.ptr(), g0_cpu_.G0_coeff_.values(),
             G0_coeff_.size() * sizeof(decltype(G0_coeff_.ptr())), cudaMemcpyHostToDevice);
  // Copy pointer to the data structure.
  BaseClass::values_ = G0_coeff_.ptr();
  BaseClass::g0_minus_ = g0_minus_.ptr();

  GlobalMemoryManager::initializeInterpolation(getStride());
  cudaDeviceSynchronize();
}

template <class InputDmn>
G0Interpolation<linalg::GPU>::G0Interpolation(const func::function<double, InputDmn>& G0_pars_t) {
  initialize(G0_pars_t);
}

double G0Interpolation<linalg::GPU>::operator()(double tau, int lindex) const {
  return details::deviceInterpolationTest(*this, tau, lindex);
}

}  // dca
}  // phys
}  // solver
}  // ctint

#endif  // DCA_HAVE_CUDA
#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_G0_INTERPOLATION_GPU_HPP
