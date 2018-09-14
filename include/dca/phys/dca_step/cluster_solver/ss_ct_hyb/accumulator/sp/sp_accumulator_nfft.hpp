// Copyright (C) 2010 Philipp Werner
//
// Integrated into DCA++ by Peter Staar (taa@zurich.ibm.com) and Bart Ydens.
//
// This class organizes the measurements in the single-site hybridization QMC integration.

/*
 * The impurity self-energy can be expressed in the following form:
 * \f{eqnarray*}{
 * \Sigma_{ab}(i \nu) = \frac12 \sum_{ij} G^{-1}_{ai}(i \nu)(U_{jb} + U_{bj}) F_{ib}^j(i \nu)
 * \f}
 * Where the impurity Green function \f$G_{ab}(\tau)\f$ is measured according to
 * \f{eqnarray*}{
 * G_{ab}(\tau) = -\frac{1}{\beta} \left< \sum_{\alpha \beta = 1}^{n} M_{\alpha \beta}
 * \delta^{-}(\tau - (\tau'_{\alpha} - \tau_{\beta}) \delta_{a,\alpha} \delta_{b,\beta} \right>_{MC}
 * \f}
 * and the correlation function \f$F_{ab}^j(\tau)\f$ is measured according to
 * \f{eqnarray*}{
 * F_{ab}^j(\tau) = -\frac{1}{\beta} \left< \sum_{\alpha \beta = 1}^{n} M_{\alpha \beta}
 * \delta^{-}(\tau - (\tau'_{\alpha} - \tau_{\beta}) n_j(\tau_{\beta}) \delta_{a,\alpha}
 * \delta_{b,\beta} \right>_{MC}
 * \f}
 * where
 * \f{eqnarray*}{
 * \delta^{\pm}(\tau) = \sum_{n} (\pm1)^n \delta(\tau - n\beta)
 * \f}
 */

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_ACCUMULATOR_SP_SP_ACCUMULATOR_NFFT_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_ACCUMULATOR_SP_SP_ACCUMULATOR_NFFT_HPP

#include <complex>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/math/nfft/dnfft_1d.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace cthyb {
// dca::phys::solver::cthyb::

template <class parameters_type, class base_cluster_type>
class SpAccumulatorNfft {
public:
  using concurrency_type = typename parameters_type::concurrency_type;
  using scalar_type = double;

  using w = func::dmn_0<domains::frequency_domain>;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using CDA = ClusterDomainAliases<parameters_type::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using KClusterDmn = typename CDA::KClusterDmn;
  // Does this _t notation carry some information?
  using r_dmn_t = RClusterDmn;
  using k_dmn_t = KClusterDmn;

  using p_dmn_t = func::dmn_variadic<nu, nu, r_dmn_t>;

  void sumTo(SpAccumulatorNfft<parameters_type, base_cluster_type>& other) const;

public:
  SpAccumulatorNfft(parameters_type& parameters_ref);

  void initialize(func::function<std::complex<double>, func::dmn_variadic<nu, nu, r_dmn_t, w>>& G_r_w,
                  func::function<std::complex<double>, func::dmn_variadic<nu, nu, r_dmn_t, w>>& GS_r_w);

  template <class configuration_type, class M_matrices_type, class H_type>
  void accumulate(double current_sign, configuration_type& configuration,
                  M_matrices_type& M_matrices, H_type& H_interactions);

  void finalize(func::function<std::complex<double>, func::dmn_variadic<nu, nu, r_dmn_t, w>>& G_r_w,
                func::function<std::complex<double>, func::dmn_variadic<nu, nu, r_dmn_t, w>>& GS_r_w);

private:
  template <class configuration_type, class H_type>
  double compute_U_times_n_2(configuration_type& configuration, H_type& H_interactions,
                             double t_start, int flavor);

private:
  parameters_type& parameters;
  concurrency_type& concurrency;

  int N_spin_orbitals;

  math::nfft::Dnfft1D<double, w, p_dmn_t> cached_nfft_1D_G_obj;
  // math::nfft::Dnfft1D<double, w, p_dmn_t> cached_nfft_1D_G_squared_obj;

  math::nfft::Dnfft1D<double, w, p_dmn_t> cached_nfft_1D_GS_obj;
  // math::nfft::Dnfft1D<double, w, p_dmn_t> cached_nfft_1D_GS_squared_obj;
};

template <class parameters_type, class base_cluster_type>
SpAccumulatorNfft<parameters_type, base_cluster_type>::SpAccumulatorNfft(parameters_type& parameters_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      N_spin_orbitals(b::dmn_size() * s::dmn_size()),

      cached_nfft_1D_G_obj(),
      cached_nfft_1D_GS_obj() {}

template <class parameters_type, class base_cluster_type>
void SpAccumulatorNfft<parameters_type, base_cluster_type>::initialize(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, r_dmn_t, w>>& G_r_w,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, r_dmn_t, w>>& GS_r_w) {
  {
    cached_nfft_1D_G_obj.resetAccumulation();
    cached_nfft_1D_GS_obj.resetAccumulation();

    G_r_w = 0;
    GS_r_w = 0;
  }
}

template <class parameters_type, class base_cluster_type>
void SpAccumulatorNfft<parameters_type, base_cluster_type>::finalize(
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, r_dmn_t, w>>& G_r_w,
    func::function<std::complex<double>, func::dmn_variadic<nu, nu, r_dmn_t, w>>& GS_r_w) {
  double beta = parameters.get_beta();

  {
    func::function<std::complex<double>, func::dmn_variadic<w, p_dmn_t>> tmp("tmp G");

    cached_nfft_1D_G_obj.finalize(tmp);

    for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++)
      for (int r_ind = 0; r_ind < r_dmn_t::dmn_size(); r_ind++)
        for (int s2_ind = 0; s2_ind < s::dmn_size(); s2_ind++)
          for (int b2_ind = 0; b2_ind < b::dmn_size(); b2_ind++)
            for (int s1_ind = 0; s1_ind < s::dmn_size(); s1_ind++)
              for (int b1_ind = 0; b1_ind < b::dmn_size(); b1_ind++)
                G_r_w(b1_ind, s1_ind, b2_ind, s2_ind, r_ind, w_ind) =
                    tmp(w_ind, b1_ind, s1_ind, b2_ind, s2_ind, r_ind);

    double one_div_n_sites = 1. / double(-beta * RClusterDmn::dmn_size());
    G_r_w *= one_div_n_sites;
  }

  {
    func::function<std::complex<double>, func::dmn_variadic<w, p_dmn_t>> tmp("tmp GS");

    cached_nfft_1D_GS_obj.finalize(tmp);

    for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++)
      for (int r_ind = 0; r_ind < r_dmn_t::dmn_size(); r_ind++)
        for (int s2_ind = 0; s2_ind < s::dmn_size(); s2_ind++)
          for (int b2_ind = 0; b2_ind < b::dmn_size(); b2_ind++)
            for (int s1_ind = 0; s1_ind < s::dmn_size(); s1_ind++)
              for (int b1_ind = 0; b1_ind < b::dmn_size(); b1_ind++)
                GS_r_w(b1_ind, s1_ind, b2_ind, s2_ind, r_ind, w_ind) =
                    tmp(w_ind, b1_ind, s1_ind, b2_ind, s2_ind, r_ind);

    double one_div_n_sites = 1. / double(-beta * RClusterDmn::dmn_size());
    GS_r_w *= one_div_n_sites;
  }
}

template <class parameters_type, class base_cluster_type>
template <class configuration_type, class M_matrices_type, class H_type>
void SpAccumulatorNfft<parameters_type, base_cluster_type>::accumulate(
    double current_sign, configuration_type& configuration, M_matrices_type& M_matrices,
    H_type& H_interactions) {
  typedef typename configuration_type::orbital_configuration_type orbital_configuration_type;

  int coor[2];
  int coor_nfft[5];

  nu nu_obj;

  scalar_type BETA = parameters.get_beta();
  scalar_type one_div_two_beta = 1. / (2. * BETA);

  double U_times_n;
  double t_i, t_j;

  for (int flavor = 0; flavor < N_spin_orbitals; flavor++) {
    nu_obj.linind_2_subind(flavor, coor);

    coor_nfft[0] = coor[0];
    coor_nfft[1] = coor[1];
    coor_nfft[2] = coor[0];
    coor_nfft[3] = coor[1];
    coor_nfft[4] = 0;

    dca::linalg::Matrix<double, dca::linalg::CPU>& M_ind = M_matrices(flavor);
    orbital_configuration_type& segments = configuration.get_vertices(flavor);

    for (int j = 0; j < (int)segments.size(); j++) {
      t_j = segments[j].t_start();
      U_times_n = compute_U_times_n_2(configuration, H_interactions, t_j, flavor);

      for (int i = 0; i < (int)segments.size(); i++) {
        t_i = segments[i].t_end();

        double scaled_tau = (t_i - t_j) * one_div_two_beta;

        cached_nfft_1D_G_obj.accumulate(coor_nfft, scaled_tau, M_ind(j, i) * current_sign);
        cached_nfft_1D_GS_obj.accumulate(coor_nfft, scaled_tau,
                                         U_times_n * M_ind(j, i) * current_sign);
      }
    }
  }
}

template <class parameters_type, class base_cluster_type>
template <class configuration_type, class H_type>
double SpAccumulatorNfft<parameters_type, base_cluster_type>::compute_U_times_n_2(
    configuration_type& configuration, H_type& H_interactions, double t_start, int flavor) {
  typedef typename configuration_type::orbital_configuration_type orbital_configuration_type;

  double U_times_n = 0;

  for (int j = 0; j < N_spin_orbitals; j++) {
    if (j == flavor)
      continue;

    if (configuration.get_full_line(j)) {
      U_times_n += 1. / 2. * (H_interactions(j, flavor, 0) + H_interactions(flavor, j, 0));
    }
    else {
      orbital_configuration_type& vertices_j = configuration.get_vertices(j);

      for (int y = 0; y < int(vertices_j.size()); y++) {
        if (vertices_j[y].t_start() < vertices_j[y].t_end()) {
          if (vertices_j[y].t_start() <= t_start && t_start <= vertices_j[y].t_end())
            U_times_n += 1. / 2. * (H_interactions(j, flavor, 0) + H_interactions(flavor, j, 0));
        }
        else {
          if ((vertices_j[y].t_start() <= t_start || t_start <= vertices_j[y].t_end()))
            U_times_n += 1. / 2. * (H_interactions(j, flavor, 0) + H_interactions(flavor, j, 0));
        }
      }
    }
  }

  return U_times_n;
}

template <class parameters_type, class base_cluster_type>
void SpAccumulatorNfft<parameters_type, base_cluster_type>::sumTo(
    SpAccumulatorNfft<parameters_type, base_cluster_type>& other) const {
  other.cached_nfft_1D_G_obj += cached_nfft_1D_G_obj;
  other.cached_nfft_1D_GS_obj += cached_nfft_1D_GS_obj;
}

}  // cthyb
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_ACCUMULATOR_SP_SP_ACCUMULATOR_NFFT_HPP
