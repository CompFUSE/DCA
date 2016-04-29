//-*-C++-*-

#ifndef SS_HYBRIDIZATION_ACCUMULATOR_SP_NFFT_H
#define SS_HYBRIDIZATION_ACCUMULATOR_SP_NFFT_H
#include "phys_library/domain_types.hpp"
#include "math_library/NFFT/dnfft_1D.h"
using namespace types;

namespace DCA {
namespace QMCI {
/*!
 *  \brief   This class organizes the measurements in the single-site hybridization QMC
 *  \author  Peter Staar
 *  \author  Bart Ydens
 *  \version 1.0
 *
 *
 *
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
 *
 */
template <class parameters_type, class base_cluster_type>
class MC_single_particle_accumulator<SS_CT_HYB, NFFT, parameters_type, base_cluster_type> {
  typedef r_DCA r_dmn_t;
  typedef k_DCA k_dmn_t;

  typedef w w_dmn_t;
  typedef dmn_3<nu, nu, r_dmn_t> p_dmn_t;

  typedef b b_dmn_t;
  typedef s s_dmn_t;

  typedef typename parameters_type::profiler_type profiler_type;
  typedef typename parameters_type::concurrency_type concurrency_type;

  typedef double scalar_type;

public:
  MC_single_particle_accumulator(parameters_type& parameters_ref);

  ~MC_single_particle_accumulator();

  void initialize(FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w>>& G_r_w,
                  FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w>>& GS_r_w);

  template <class walker_type, class H_type>
  void accumulate(walker_type& walker, H_type& H_interactions);

  template <class configuration_type, class M_matrices_type, class H_type>
  void accumulate(double current_sign, configuration_type& configuration,
                  M_matrices_type& M_matrices, H_type& H_interactions);

  void finalize(FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w>>& G_r_w,
                FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w>>& GS_r_w);

private:
  template <class walker_type, class H_type>
  double compute_U_times_n(walker_type& walker, H_type& H_interactions, double t_start, int flavor);

  template <class configuration_type, class H_type>
  double compute_U_times_n_2(configuration_type& configuration, H_type& H_interactions,
                             double t_start, int flavor);

private:
  parameters_type& parameters;
  concurrency_type& concurrency;

  int N_spin_orbitals;

  math_algorithms::NFFT::dnfft_1D<double, w_dmn_t, p_dmn_t> cached_nfft_1D_G_obj;
  // math_algorithms::NFFT::dnfft_1D<double, w_dmn_t, p_dmn_t> cached_nfft_1D_G_squared_obj;

  math_algorithms::NFFT::dnfft_1D<double, w_dmn_t, p_dmn_t> cached_nfft_1D_GS_obj;
  // math_algorithms::NFFT::dnfft_1D<double, w_dmn_t, p_dmn_t> cached_nfft_1D_GS_squared_obj;
};

template <class parameters_type, class base_cluster_type>
MC_single_particle_accumulator<SS_CT_HYB, NFFT, parameters_type, base_cluster_type>::MC_single_particle_accumulator(
    parameters_type& parameters_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      N_spin_orbitals(b::dmn_size() * s::dmn_size()),

      cached_nfft_1D_G_obj(),
      cached_nfft_1D_GS_obj() {}

template <class parameters_type, class base_cluster_type>
MC_single_particle_accumulator<SS_CT_HYB, NFFT, parameters_type,
                               base_cluster_type>::~MC_single_particle_accumulator() {}

template <class parameters_type, class base_cluster_type>
void MC_single_particle_accumulator<SS_CT_HYB, NFFT, parameters_type, base_cluster_type>::initialize(
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w>>& G_r_w,
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w>>& GS_r_w) {
  {
    cached_nfft_1D_G_obj.initialize();
    cached_nfft_1D_GS_obj.initialize();

    G_r_w = 0;
    GS_r_w = 0;
  }
}

template <class parameters_type, class base_cluster_type>
void MC_single_particle_accumulator<SS_CT_HYB, NFFT, parameters_type, base_cluster_type>::finalize(
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w>>& G_r_w,
    FUNC_LIB::function<std::complex<double>, dmn_4<nu, nu, r_dmn_t, w>>& GS_r_w) {
  double beta = parameters.get_beta();

  {
    FUNC_LIB::function<std::complex<double>, dmn_2<w_dmn_t, p_dmn_t>> tmp("tmp G");

    cached_nfft_1D_G_obj.finalize(tmp);

    for (int w_ind = 0; w_ind < w_dmn_t::dmn_size(); w_ind++)
      for (int r_ind = 0; r_ind < r_dmn_t::dmn_size(); r_ind++)
        for (int s2_ind = 0; s2_ind < s_dmn_t::dmn_size(); s2_ind++)
          for (int b2_ind = 0; b2_ind < b_dmn_t::dmn_size(); b2_ind++)
            for (int s1_ind = 0; s1_ind < s_dmn_t::dmn_size(); s1_ind++)
              for (int b1_ind = 0; b1_ind < b_dmn_t::dmn_size(); b1_ind++)
                G_r_w(b1_ind, s1_ind, b2_ind, s2_ind, r_ind, w_ind) =
                    tmp(w_ind, b1_ind, s1_ind, b2_ind, s2_ind, r_ind);

    double one_div_n_sites = 1. / double(-beta * r_DCA::dmn_size());
    G_r_w *= one_div_n_sites;
  }

  {
    FUNC_LIB::function<std::complex<double>, dmn_2<w_dmn_t, p_dmn_t>> tmp("tmp GS");

    cached_nfft_1D_GS_obj.finalize(tmp);

    for (int w_ind = 0; w_ind < w_dmn_t::dmn_size(); w_ind++)
      for (int r_ind = 0; r_ind < r_dmn_t::dmn_size(); r_ind++)
        for (int s2_ind = 0; s2_ind < s_dmn_t::dmn_size(); s2_ind++)
          for (int b2_ind = 0; b2_ind < b_dmn_t::dmn_size(); b2_ind++)
            for (int s1_ind = 0; s1_ind < s_dmn_t::dmn_size(); s1_ind++)
              for (int b1_ind = 0; b1_ind < b_dmn_t::dmn_size(); b1_ind++)
                GS_r_w(b1_ind, s1_ind, b2_ind, s2_ind, r_ind, w_ind) =
                    tmp(w_ind, b1_ind, s1_ind, b2_ind, s2_ind, r_ind);

    double one_div_n_sites = 1. / double(-beta * r_DCA::dmn_size());
    GS_r_w *= one_div_n_sites;
  }
}

template <class parameters_type, class base_cluster_type>
template <class configuration_type, class M_matrices_type, class H_type>
void MC_single_particle_accumulator<SS_CT_HYB, NFFT, parameters_type, base_cluster_type>::accumulate(
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

    LIN_ALG::matrix<double, LIN_ALG::CPU>& M_ind = M_matrices(flavor);
    orbital_configuration_type& segments = configuration.get_vertices(flavor);

    for (int j = 0; j < (int)segments.size(); j++) {
      t_j = segments[j].t_start();
      U_times_n = compute_U_times_n_2(configuration, H_interactions, t_j, flavor);

      for (int i = 0; i < (int)segments.size(); i++) {
        t_i = segments[i].t_end();

        double scaled_tau = (t_i - t_j) * one_div_two_beta;

        cached_nfft_1D_G_obj.accumulate_at(coor_nfft, scaled_tau, M_ind(j, i) * current_sign);
        cached_nfft_1D_GS_obj.accumulate_at(coor_nfft, scaled_tau,
                                            U_times_n * M_ind(j, i) * current_sign);
      }
    }
  }
}

template <class parameters_type, class base_cluster_type>
template <class configuration_type, class H_type>
double MC_single_particle_accumulator<SS_CT_HYB, NFFT, parameters_type, base_cluster_type>::compute_U_times_n_2(
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
}
}

#endif
