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

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_G0_INTERPOLATION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_G0_INTERPOLATION_HPP

#include <assert.h>
#include <vector>

#include "dca/function/domains/dmn.hpp"
#include "dca/function/domains/dmn_0.hpp"
#include "dca/function/domains/dmn_variadic.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/device_type.hpp"
#include "dca/math/interpolation/akima_interpolation.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/domains/ct_int_domains.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/function_proxy.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

// void template
template <linalg::DeviceType device_t>
class G0Interpolation {};

// ParametersDomain is a collection of discrete labels not involving the time.
template <>
class G0Interpolation<linalg::CPU> {
private:
  using Pdmn = G0ParametersDomain;
  using Pdmn0 = func::dmn_0<G0ParametersDomain>;
  using Tdmn = func::dmn_0<domains::time_domain>;
  using PTdmn = func::dmn_variadic<Pdmn0, Tdmn>;

public:
  G0Interpolation() = default;
  // See this->initialize(G0_pars_t).
  template <class InputDmn>
  G0Interpolation(const func::function<double, InputDmn>& G0_pars_t);

  virtual ~G0Interpolation() {}

  // In: G0_pars_t. Assumed to be a function of discrete labels and time (in this order) and to
  //             be antiperiodic in time.
  template <class InputDmn>
  void initialize(const func::function<double, InputDmn>& G0_pars_t);
  // Constructor. Reshape the G0Dmn and calls the first constructor.
  // In: G0_r_t: function of dmn_variadic<D1, D2, ..., Tdmn>

  // Returns cubic interpolation of G0(tau) in the spin-band-position defined by lindex.
  double operator()(double tau, int lindex) const;

  // Number of value if g0 stored per parameter value.
  int getTimeSlices() const {
    return G0_coeff_[1];
  }
  int getStride() const {
    return getTimeSlices() * COEFF_SIZE;
  }

  friend class G0Interpolation<linalg::GPU>;

  static constexpr int COEFF_SIZE = 4;

private:
  virtual void initialize(const FunctionProxy<double, PTdmn>& G0_pars_t);

private:
  using CoeffDmn = func::dmn_0<func::dmn<COEFF_SIZE>>;
  using PTime0 = func::dmn_0<PositiveTimeDomain>;
  using InterpolationDmn = func::dmn_variadic<CoeffDmn, PTime0, Pdmn0>;

  func::function<double, InterpolationDmn> G0_coeff_;
  double beta_ = 0;
  // value at tau = 0
  std::vector<double> g0_minus_;
  // Spacing between time bins.
  double n_div_beta_;
};

template <class InputDmn>
G0Interpolation<linalg::CPU>::G0Interpolation(const func::function<double, InputDmn>& G0_pars_t) {
  initialize(G0_pars_t);
}

template <class InputDmn>
void G0Interpolation<linalg::CPU>::initialize(const func::function<double, InputDmn>& G0_pars_t) {
  PositiveTimeDomain::initialize();
  Pdmn::initialize(InputDmn::dmn_size() / Tdmn::dmn_size());
  initialize(FunctionProxy<double, PTdmn>(G0_pars_t));
}

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTINT_WALKER_TOOLS_G0_INTERPOLATION_HPP
