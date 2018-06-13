// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
// See LICENSE for terms of usage./
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides type definitions and helper methods for the ED cluster solver tests.

#ifndef TEST_INTEGRATION_EXACT_DIAGONALIZATION_ADVANCED_ED_CLUSTER_SOLVER_TEST_HELPER_HPP
#define TEST_INTEGRATION_EXACT_DIAGONALIZATION_ADVANCED_ED_CLUSTER_SOLVER_TEST_HELPER_HPP

#include <cassert>
#include <cmath>
#include <complex>

#include "dca/function/function.hpp"
#include "dca/function/domains.hpp"
#include "dca/phys/dca_data/dca_data.hpp"
#include "dca/phys/dca_data/dca_data_real_freq.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include "dca/phys/models/analytic_hamiltonians/square_lattice.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/parameters.hpp"

namespace dca {
namespace testing {
// dca::testing::

// Typedefs for domains
using MatsubaraFreqDmn = func::dmn_0<phys::domains::frequency_domain>;
using ImagTimeDmn = func::dmn_0<phys::domains::time_domain>;
using OrbitalSpinDmn = func::dmn_variadic<func::dmn_0<phys::domains::electron_band_domain>,
                                          func::dmn_0<phys::domains::electron_spin_domain>>;
using KDmn = func::dmn_0<phys::domains::cluster_domain<
    double, 2, phys::domains::CLUSTER, phys::domains::MOMENTUM_SPACE, phys::domains::BRILLOUIN_ZONE>>;
using RDmn = func::dmn_0<phys::domains::cluster_domain<
    double, 2, phys::domains::CLUSTER, phys::domains::REAL_SPACE, phys::domains::BRILLOUIN_ZONE>>;

// Typedefs for common dca::func::functions
using F_k_w =
    func::function<std::complex<double>,
                   func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KDmn, MatsubaraFreqDmn>>;
using F_k_t =
    func::function<double, func::dmn_variadic<OrbitalSpinDmn, OrbitalSpinDmn, KDmn, ImagTimeDmn>>;

// Typedefs for model, parameters and data
using Lattice = phys::models::square_lattice<phys::domains::D4>;
using Model = phys::models::TightBindingModel<Lattice>;
using Parameters = phys::params::Parameters<parallel::NoConcurrency, void, void, Model, void,
                                            phys::solver::CT_AUX>;  // CT_AUX is a placeholder
using Data = phys::DcaData<Parameters>;
using DataRealFreq = phys::DcaDataRealFreq<Parameters>;

// Computes the analytic form of the free Matsubara Green's function,
// G_0(\vec{k}, i\omega) = \left( i\omega - \epsilon_{\vec{k}} + \mu \right)^{-1}.
void computeAnalyticG0_k_w(
    /*const*/ func::function<double, func::dmn_variadic<OrbitalSpinDmn, KDmn>>& eps_k,
    const double mu, F_k_w& G0_k_w) {
  const std::complex<double> i(0., 1.);

  for (int w = 0; w < MatsubaraFreqDmn::dmn_size(); ++w) {
    const double w_val = MatsubaraFreqDmn::get_elements()[w];
    for (int k = 0; k < KDmn::dmn_size(); ++k) {
      for (int m = 0; m < OrbitalSpinDmn::dmn_size(); ++m) {
        G0_k_w(m, m, k, w) = 1. / (i * w_val - eps_k(m, k) + mu);
      }
    }
  }
}

// Computes the analytic form of the free imaginary time Green's function,
// G_0(\vec{k}, \tau >= 0) =
//     - \frac{e^{-\tau (\epsilon_{\vec{k}}-\mu)}}{e^{-\beta(\epsilon_{\vec{k}}-\mu)}+1},
//         if (\epsilon_{\vec{k}}-\mu) >= 0,
//    - \frac{e^{(\beta-\tau) (\epsilon_{\vec{k}}-\mu)}}{e^{\beta(\epsilon_{\vec{k}}-\mu)}+1},
//        if (\epsilon_{\vec{k}}-\mu) < 0,
// G_0(\vec{k}, \tau < 0) = -G_0(\vec{k}, \tau + \beta).
void computeAnalyticG0_k_t(
    /*const*/ func::function<double, func::dmn_variadic<OrbitalSpinDmn, KDmn>>& eps_k,
    const double mu, const double beta, F_k_t& G0_k_t) {
  const int t0 = ImagTimeDmn::dmn_size() / 2;

  for (int t = 0; t < t0; ++t) {
    const double tau = ImagTimeDmn::get_elements()[t + t0];
    for (int k = 0; k < KDmn::dmn_size(); ++k) {
      for (int m = 0; m < OrbitalSpinDmn::dmn_size(); ++m) {
        const double eps_min_mu = eps_k(m, k) - mu;
        const double val =
            eps_min_mu >= 0.
                ? -std::exp(-tau * eps_min_mu) / (std::exp(-beta * eps_min_mu) + 1)
                : -std::exp((beta - tau) * eps_min_mu) / (std::exp(beta * eps_min_mu) + 1);
        G0_k_t(m, m, k, t + t0) = val;
        G0_k_t(m, m, k, t) = -val;
      }
    }
  }
}

// Computes the analytic form of the Matsubara Green's function of a single-site Hubbard model at
// half-filling,
// G(i\omega) = \left( i \omega - \frac{U^2}{4 i\omega} \right)^{-1}.
void computeAnalyticG_w(const double U, F_k_w& G_w) {
  const std::complex<double> i(0., 1.);

  for (int m = 0; m < OrbitalSpinDmn::dmn_size(); ++m) {
    for (int w = 0; w < MatsubaraFreqDmn::dmn_size(); ++w) {
      const double w_val = MatsubaraFreqDmn::get_elements()[w];
      G_w(m, m, 0, w) = 1. / (i * w_val - U * U / (4. * i * w_val));
    }
  }
}

}  // testing
}  // dca

#endif  // TEST_INTEGRATION_EXACT_DIAGONALIZATION_ADVANCED_ED_CLUSTER_SOLVER_TEST_HELPER_HPP
