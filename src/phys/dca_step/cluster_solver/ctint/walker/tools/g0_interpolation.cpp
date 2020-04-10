// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
//  See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// CPU implementation of g0_interpolation.hpp.

#include <dca/phys/dca_step/cluster_solver/ctint/walker/tools/g0_interpolation.hpp>

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

template <typename Real>
void G0Interpolation<linalg::CPU, Real>::initialize(const FunctionProxy<double, PTdmn>& G0_pars_t) {
  beta_ = PositiveTimeDomain::get_elements().back();
  n_div_beta_ = Real(PositiveTimeDomain::get_size() - 1) / beta_;

  const int t_pos_size = PositiveTimeDomain::get_size();
  dca::math::interpolation::akima_interpolation<Real> akima_obj(t_pos_size);
  G0_coeff_.reset();
  g0_minus_.resize(Pdmn::get_size());

  std::vector<Real> y(t_pos_size), x(t_pos_size);
  // The abscissa is scaled to integer steps
  for (int i = 0; i < x.size(); i++)
    x[i] = i;

  // Loop over each cluster-orbital label.
  for (int p = 0; p < Pdmn::get_size(); p++) {
    // set G0(0-)
    g0_minus_[p] = G0_pars_t(p, t_pos_size - 1);
    // Compute interpolation coefficients:
    for (int t = 0; t < t_pos_size; t++)
      y[t] = G0_pars_t(p, t + t_pos_size);
    // INTERNAL initialize or initialize_periodic ?
    akima_obj.initialize(x.data(), y.data());
    // Store coefficients:
    for (int t = 0; t < t_pos_size - 1; t++) {
      for (int l = 0; l < 4; l++)
        G0_coeff_(l, t, p) = akima_obj.get_alpha(l, t);
      assert(std::abs(G0_coeff_(0, t, p) - G0_pars_t(p, t + t_pos_size)) <
             10 * std::numeric_limits<Real>::epsilon());
    }
  }
}

template <typename Real>
Real G0Interpolation<linalg::CPU, Real>::operator()(Real tau, int lindex) const {
  assert(beta_ != 0);
  if (tau == 0)  // returns G0(tau = 0+)
    return g0_minus_[lindex];

  short int factor = 1;
  if (tau < 0) {
    tau += beta_;
    factor = -1;
  }
  assert(tau >= 0 and tau <= beta_);
  assert(lindex >= 0 and lindex < Pdmn::get_size());

  // Scale tau in [0, n_time_slices). Assume even spacing in time.
  const Real scaled_tau = tau * n_div_beta_;
  const int tau_index(scaled_tau);
  const Real delta_tau = scaled_tau - tau_index;

  // Get the pointer to the first akima coeff.
  const Real* const coeff_ptr = &G0_coeff_(0, tau_index, lindex);
  // Return akima interpolation.
  return factor *
         (coeff_ptr[0] +
          delta_tau * (coeff_ptr[1] + delta_tau * (coeff_ptr[2] + delta_tau * coeff_ptr[3])));
}

// Instantation
template class G0Interpolation<linalg::CPU, float>;
template class G0Interpolation<linalg::CPU, double>;

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca
