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

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_INTERPOLATION_G0_INTERPOLATION_CPU_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_INTERPOLATION_G0_INTERPOLATION_CPU_HPP

#include <assert.h>
#include <vector>

#include "dca/function/domains/dmn.hpp"
#include "dca/function/domains/dmn_0.hpp"
#include "dca/function/domains/dmn_variadic.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/device_type.hpp"
#include "dca/math/interpolation/akima_interpolation.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/interpolation/interpolation_domains.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/interpolation/function_proxy.hpp"
#include "dca/phys/dca_step/cluster_solver/shared_tools/interpolation/shrink_G0.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include "dca/util/type_utils.hpp"

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

// void template
template <linalg::DeviceType device_t, typename Scalar>
class G0Interpolation {};

// ParametersDomain is a collection of discrete labels not involving the time.
template <typename Scalar>
class G0Interpolation<linalg::CPU, Scalar> {
protected:
  using Pdmn = G0ParametersDomain;
  using Pdmn0 = func::dmn_0<G0ParametersDomain>;
  using Tdmn = func::dmn_0<domains::time_domain>;
  using PTdmn = func::dmn_variadic<Pdmn0, Tdmn>;
  using Real = dca::util::Real<Scalar>;

public:
  G0Interpolation() = default;

  template <class Dmn, typename Scalar2>
  G0Interpolation(const func::function<Scalar2, Dmn>& G0_func) {
    initialize(G0_func);
  }

  virtual ~G0Interpolation() = default;

  // In: G0_pars_t. Assumed to be a function of discrete labels and time (in this order) and to
  //             be antiperiodic in time.
  template <class Dmn>
  void initialize(const func::function<Scalar, Dmn>& G0_func);

  template <class Dmn, typename Scalar2>
  void initialize(const func::function<Scalar2, Dmn>& G0_func) {
    const func::function<Scalar, Dmn> cpy(G0_func);
    initialize(cpy);
  }

  // Initialize with only one spin sector.
  template <int dim, typename Scalar2>
  void initializeShrinked(const details::SpGreensFunction<dim, Scalar2>& g0_r_t) {
    initialize(details::shrinkG0(g0_r_t));
  }

  // Returns cubic interpolation of G0(tau) in the spin-band-position defined by lindex.
  Scalar operator()(Real tau, int lindex) const;

  // Number of value if g0 stored per parameter value.
  int getTimeSlices() const {
    return G0_coeff_[1];
  }
  int getStride() const {
    return getTimeSlices() * COEFF_SIZE;
  }

  friend class G0Interpolation<linalg::GPU, Scalar>;

  static constexpr int COEFF_SIZE = 4;

private:
  virtual void initialize(const FunctionProxy<Scalar, PTdmn>& G0_pars_t);

private:
  using CoeffDmn = func::dmn_0<func::dmn<COEFF_SIZE>>;
  using PTime0 = func::dmn_0<PositiveTimeDomain>;
  using InterpolationDmn = func::dmn_variadic<CoeffDmn, PTime0, Pdmn0>;

  func::function<Scalar, InterpolationDmn> G0_coeff_;
  Real beta_ = 0;
  // value at tau = 0
  std::vector<Scalar> g0_minus_;
  // Spacing between time bins.
  Real n_div_beta_;
};

template <typename Scalar>
template <class Dmn>
void G0Interpolation<linalg::CPU, Scalar>::initialize(const func::function<Scalar, Dmn>& G0_func) {
  PositiveTimeDomain::initialize();
  Pdmn::initialize(Dmn::dmn_size() / Tdmn::dmn_size());
  initialize(FunctionProxy<Scalar, PTdmn>(G0_func));
}

}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SHARED_TOOLS_INTERPOLATION_G0_INTERPOLATION_CPU_HPP
