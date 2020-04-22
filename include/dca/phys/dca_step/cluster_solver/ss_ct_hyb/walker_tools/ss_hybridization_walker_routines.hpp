// Copyright (C) 2010 Philipp Werner
//
// Integrated into DCA++ by Peter Staar (taa@zurich.ibm.com) and Bart Ydens.
//
// This class implements the helper functions for the insertion and removal of (anti-)segments. The
// helper functions include the calculation of the determinant ratio and the computation of the new
// hybridization matrix using sherman-morrison equations.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_WALKER_TOOLS_SS_HYBRIDIZATION_WALKER_ROUTINES_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_WALKER_TOOLS_SS_HYBRIDIZATION_WALKER_ROUTINES_HPP

#include <cassert>
#include <iostream>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/linalg.hpp"
#include "dca/math/interpolation/akima_interpolation.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/ss_hybridization_solver_routines.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/structures/hybridization_vertex.hpp"
#include "dca/phys/dca_step/cluster_solver/ss_ct_hyb/walker_tools/static_matrix_routines.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain_left_oriented.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace cthyb {
// dca::phys::solver::cthyb::

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
class ss_hybridization_walker_routines
    : public ss_hybridization_solver_routines<parameters_t, MOMS_t> {
public:
  using t = func::dmn_0<domains::time_domain>;
  using w = func::dmn_0<domains::frequency_domain>;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  typedef parameters_t parameters_type;
  using CDA = ClusterDomainAliases<parameters_type::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;

  using nu_nu_r_DCA_t = func::dmn_variadic<nu, nu, RClusterDmn, t>;

  typedef func::dmn_0<domains::time_domain_left_oriented> shifted_t;
  typedef func::dmn_variadic<nu, nu, RClusterDmn, shifted_t> nu_nu_r_dmn_t_shifted_t;

  typedef func::dmn_0<func::dmn<4, int>> akima_dmn_t;
  typedef func::dmn_variadic<akima_dmn_t, nu, nu, RClusterDmn, shifted_t> akima_nu_nu_r_dmn_t_shifted_t;

  typedef MOMS_t MOMS_type;
  typedef configuration_t configuration_type;
  typedef rng_t rng_type;

  typedef typename parameters_type::concurrency_type concurrency_type;

public:
  ss_hybridization_walker_routines(parameters_t& parameters_ref, MOMS_t& MOMS_ref,
                                   configuration_t& configuration_ref, rng_t& rng_ref);

  void initialize() {}
  void initialize_akima_coefficients(func::function<double, nu_nu_r_DCA_t>& F_r_t);

  const parameters_t& get_parameters() {
    return parameters;
  }
  MOMS_t& get_MOMS() {
    return MOMS;
  }
  configuration_t& get_configuration() {
    return configuration;
  }
  rng_t& get_rng() {
    return rng;
  }

  static int cycle(int i, int size);

  // matrix-routines --> has to move to LIN_ALG at some point!
  template <typename scalar_type>
  static void cycle_column_forward(dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& M);

  template <typename scalar_type>
  static void cycle_column_backward(dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& M);

  template <typename scalar_type>
  static void cycle_row_forward(dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& M);

  template <typename scalar_type>
  static void cycle_row_backward(dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& M);

  template <typename Hybridization_function_t>
  double interpolate_F(int* coor_flavor, double tau, Hybridization_function_t& F);

  double compute_length(double r, double l_max, double mu);

  template <typename S>
  double compute_overlap(Hybridization_vertex segment, S& other_segments, int other_full_line,
                         double BETA);

  template <typename S>
  double segment_overlap(Hybridization_vertex segment, S& other_segments, int other_full_line,
                         double BETA);

  template <typename orbital_configuration_t>
  void compute_intervals(double t, double BETA, double& t_up, double& t_down,
                         orbital_configuration_t& segments,
                         typename orbital_configuration_t::iterator& s_up,
                         typename orbital_configuration_t::iterator& s_down);

  // Calculates the determinant ratio for inserting a new vertex. The determinant ratio is given by
  // (A.10)
  // \f{eqnarray*}{
  // \frac{det(M^{k+1}_{\sigma})}{det(M^{k}_{\sigma})} = det(S - R M^k Q)
  // \f}
  // with\f$ S = F(\tau^n_{e} - \tau^n_{s})\f$,
  // R a (1 x k)-vector with \f$R[i] =  F(\tau^n_{e} - \tau^i_{s})\f$,
  // Q a (k x 1)-vector with \f$Q[i] =  F(\tau^i_{e} - \tau^n_{s})\f$.
  template <typename G, typename vertex_vertex_matrix_type, typename orbital_configuration_t>
  double det_rat_up(int this_flavor, Hybridization_vertex& new_segment, vertex_vertex_matrix_type& M,
                    orbital_configuration_t& segments_old, G& F, std::vector<double>& R,
                    std::vector<double>& Q, double& det_rat_sign, double& overlap);

  // Calculates the new hybridization matrix for inserting a new vertex using sherman-morrison
  // equations (A.4-9).
  template <typename G, typename vertex_vertex_matrix_type, typename orbital_configuration_t>
  void compute_M_up(int r, int s, vertex_vertex_matrix_type& M, orbital_configuration_t& segments_old,
                    G& F, std::vector<double>& Fs, std::vector<double>& Fe, double det_rat);

  // Calculates the determinant ratio for removing a vertex.
  template <typename vertex_vertex_matrix_type, typename orbital_configuration_t>
  double det_rat_down(int r, int s, vertex_vertex_matrix_type& M,
                      orbital_configuration_t& segments_old, double& det_rat_sign);

  // Calculates the new hybridization matrix for removing a vertex using sherman-morrison equations
  // (A.4-9).
  template <typename vertex_vertex_matrix_type>
  void compute_M_down(int r, int s, vertex_vertex_matrix_type& M);

  // Calculates the determinant ratio for shifting a vertex end point. (u is k-th unity vector)
  // \f{eqnarray}{
  // det_rat &=& 1 + v*A^{-1}*u
  //         &=& 1 + (F_{new} - F_{old}) * A^{-1} *u
  //         &=& 1 + F_{new} * A^{-1} *u -  F_{old} * A^{-1} *u = F_{new} * A^{-1} *u
  // \f}
  // \f$ F_{old} \f$ is k-th row of matrix A, and \f$A^{-1} *u\f$ is k_th column of \f$A^{-1}\f$ and
  // thus \f$ F_{old} *A^{-1} *u\f$ is equal to 1. (\f$A A^{-1} = I\f$)

  template <typename G, typename vertex_vertex_matrix_type, typename orbital_configuration_t>
  double det_rat_shift_end(int this_flavor, Hybridization_vertex& new_segment, int k,
                           vertex_vertex_matrix_type& M, orbital_configuration_t& segments_old,
                           G& F, std::vector<double>& R, std::vector<double>& Q,
                           double& det_rat_sign, double& overlap);

  // Calculates the determinant ratio for shifting a vertex start point. (v is k-th unity vector)
  // \f{eqnarray}{
  // det_rat &=& 1 + v*A^{-1}*u
  //         &=& 1 + v * A^{-1} *(F_{new} - F_{old})
  //         &=& 1 + v * A^{-1} *F_{new} -  v * A^{-1} * F_{old} = v * A^{-1} *F_{new}
  // \f}
  // \f$ F_{old} \f$ is k-th column of matrix A, and \f$ v * A^{-1} \f$ is k_th row of \f$ A^{-1}
  // \f$ and thus \f$ v * A^{-1} * F_{old} \f$ is equal to 1. (\f$A A^{-1} = I\f$)
  template <typename G, typename vertex_vertex_matrix_type, typename orbital_configuration_t>
  double det_rat_shift_start(int this_flavor, Hybridization_vertex& new_segment, int k,
                             vertex_vertex_matrix_type& M, orbital_configuration_t& segments_old,
                             G& F, std::vector<double>& R, std::vector<double>& Q,
                             double& det_rat_sign, double& overlap);

  // Calculates the new hybridization matrix for shifting a vertex end point using sherman-morrison
  // equations (A.4-9). Q_prime is actually -Q_prime
  template <typename vertex_vertex_matrix_type>
  void compute_M_shift_end(int k, vertex_vertex_matrix_type& M, std::vector<double>& Fs,
                           std::vector<double>& Fe, double det_rat);

  // Calculates the new hybridization matrix for shifting a vertex start point using
  // sherman-morrison equations (A.4-9). R_prime is actually -R_prime
  template <typename vertex_vertex_matrix_type>
  void compute_M_shift_start(int k, vertex_vertex_matrix_type& M, std::vector<double>& Fs,
                             std::vector<double>& Fe, double det_rat);
  // Calculates Q' = -M*Q.
  template <typename vertex_vertex_matrix_type>
  void compute_Q_prime(std::vector<double>& Q, vertex_vertex_matrix_type& M,
                       std::vector<double>& Q_prime);

  // Calculates R' = -R*M.
  template <typename vertex_vertex_matrix_type>
  void compute_R_prime(std::vector<double>& R, vertex_vertex_matrix_type& M,
                       std::vector<double>& R_prime);

  // Calculates new M matrix M_new = Q'*R'*S'.
  template <typename vertex_vertex_matrix_type>
  void compute_M(std::vector<double>& Q_prime, std::vector<double>& R_prime, double S_prime,
                 vertex_vertex_matrix_type& M);

private:
  const parameters_t& parameters;
  MOMS_t& MOMS;
  const concurrency_type& concurrency;

  configuration_t& configuration;
  rng_t& rng;

  nu_nu_r_dmn_t_shifted_t nu_nu_r_dmn_t_t_shifted_dmn;
  func::function<double, akima_nu_nu_r_dmn_t_shifted_t> akima_coefficients;
};

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::ss_hybridization_walker_routines(
    parameters_t& parameters_ref, MOMS_t& MOMS_ref, configuration_t& configuration_ref, rng_t& rng_ref)
    : ss_hybridization_solver_routines<parameters_t, MOMS_t>(parameters_ref, MOMS_ref),

      parameters(parameters_ref),
      MOMS(MOMS_ref),
      concurrency(parameters.get_concurrency()),

      configuration(configuration_ref),
      rng(rng_ref) {
  // std::cout << __FUNCTION__ << endl;

  initialize();
}

template <class parameters_t, class MOMS_t, typename configuration_t, typename rng_t>
void ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::initialize_akima_coefficients(
    func::function<double, nu_nu_r_DCA_t>& F_r_t) {
  int size = t::dmn_size() / 2;

  math::interpolation::akima_interpolation<double> ai_obj(size);

  double* x = new double[size];
  double* y = new double[size];

  for (int t_ind = 0; t_ind < t::dmn_size() / 2; t_ind++)
    x[t_ind] = t_ind;

  {
    for (int r_ind = 0; r_ind < RClusterDmn::dmn_size(); r_ind++) {
      for (int nu1_ind = 0; nu1_ind < nu::dmn_size(); nu1_ind++) {
        for (int nu0_ind = 0; nu0_ind < nu::dmn_size(); nu0_ind++) {
          for (int t_ind = 0; t_ind < t::dmn_size() / 2; t_ind++)
            y[t_ind] = -F_r_t(nu0_ind, nu1_ind, r_ind, t_ind);

          ai_obj.initialize(x, y);

          for (int t_ind = 0; t_ind < t::dmn_size() / 2 - 1; t_ind++)
            for (int l = 0; l < 4; l++)
              akima_coefficients(l, nu0_ind, nu1_ind, r_ind, t_ind) = ai_obj.get_alpha(l, t_ind);
        }
      }
    }
  }

  {
    for (int r_ind = 0; r_ind < RClusterDmn::dmn_size(); r_ind++) {
      for (int nu1_ind = 0; nu1_ind < nu::dmn_size(); nu1_ind++) {
        for (int nu0_ind = 0; nu0_ind < nu::dmn_size(); nu0_ind++) {
          for (int t_ind = t::dmn_size() / 2; t_ind < t::dmn_size(); t_ind++)
            y[t_ind - t::dmn_size() / 2] = -F_r_t(nu0_ind, nu1_ind, r_ind, t_ind);

          ai_obj.initialize(x, y);

          for (int t_ind = t::dmn_size() / 2; t_ind < t::dmn_size() - 1; t_ind++)
            for (int l = 0; l < 4; l++)
              akima_coefficients(l, nu0_ind, nu1_ind, r_ind, t_ind - 1) =
                  ai_obj.get_alpha(l, t_ind - t::dmn_size() / 2);
        }
      }
    }
  }

  delete[] x;
  delete[] y;
}

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
template <typename Hybridization_function_t>
double ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::interpolate_F(
    int* coor_flavor, double tau, Hybridization_function_t& /*F*/) {
  const double beta = parameters.get_beta();
  const double N_div_beta = parameters.get_sp_time_intervals() / beta;

  int b_i = coor_flavor[0];
  int s_i = coor_flavor[1];
  int b_j = coor_flavor[0];
  int s_j = coor_flavor[1];
  int delta_r = 0;  // coor_flavor[4];

  // make sure that new_tau is positive !!
  // double new_tau = tau+beta;
  double new_tau = -tau + beta;
  if (new_tau < -1.e-6 and new_tau > 2 * beta + 1.e-6)
    std::cout << "\t" << tau << "\t" << new_tau << "\n";

  assert(new_tau > -1.e-6 and new_tau < 2 * beta + 1.e-6);

  double scaled_tau = new_tau * N_div_beta;

  int t_ind = scaled_tau;
  // assert(shifted_t::get_elements()[t_ind]<=tau &&
  // tau<shifted_t::get_elements()[t_ind]+1./N_div_beta);
  assert(shifted_t::get_elements()[t_ind] <= -tau &&
         -tau < shifted_t::get_elements()[t_ind] + 1. / N_div_beta);

  double delta_tau = scaled_tau - t_ind;
  assert(delta_tau > -1.e-16 && delta_tau <= 1 + 1.e-16);

  int linind = 4 * nu_nu_r_dmn_t_t_shifted_dmn(b_i, s_i, b_j, s_j, delta_r, t_ind);

  double* a_ptr = &akima_coefficients(linind);

  double result = (a_ptr[0] + delta_tau * (a_ptr[1] + delta_tau * (a_ptr[2] + delta_tau * a_ptr[3])));

  // return result;
  return -result;
}

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
int ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::cycle(int i,
                                                                                          int size) {
  return (i > 0 ? i - 1 : size - 1);
}

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
template <typename scalar_type>
void ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::cycle_column_forward(
    dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& M) {
  static_matrix_routines::cycle_column_forward(M);
}

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
template <typename scalar_type>
void ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::cycle_column_backward(
    dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& M) {
  static_matrix_routines::cycle_column_backward(M);
}

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
template <typename scalar_type>
void ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::cycle_row_forward(
    dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& M) {
  static_matrix_routines::cycle_row_forward(M);
}

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
template <typename scalar_type>
void ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::cycle_row_backward(
    dca::linalg::Matrix<scalar_type, dca::linalg::CPU>& M) {
  static_matrix_routines::cycle_row_backward(M);
}

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
double ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::compute_length(
    double r, double l_max, double mu) {
  if (mu == 0)
    return r * l_max;
  else
    return 1 / mu * log(r * (exp(mu * l_max) - 1) + 1);
}

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
template <typename S>
double ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::compute_overlap(
    Hybridization_vertex segment, S& other_segments, int other_full_line, double BETA) {
  if (segment.t_start() < segment.t_end()) {
    return segment_overlap(segment, other_segments, other_full_line, BETA);
  }
  else {
    double other_length = 0;

    Hybridization_vertex segment1(0, segment.t_end());
    Hybridization_vertex segment2(segment.t_start(), BETA);

    other_length += segment_overlap(segment1, other_segments, other_full_line, BETA);
    other_length += segment_overlap(segment2, other_segments, other_full_line, BETA);

    return other_length;
  }
}

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
template <typename S>
double ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::segment_overlap(
    Hybridization_vertex segment, S& other_segments, int other_full_line, double BETA) {
  double length = (segment.t_start() < segment.t_end() ? segment.t_end() - segment.t_start()
                                                       : segment.t_end() - segment.t_start() + BETA);
  double t_final = segment.t_start() + length;
  double t = segment.t_start();
  double t_final_segment;
  double other_length = 0;
  if (other_full_line == 1)
    other_length = length;
  else if (other_segments.size() > 0) {
    typename S::iterator it;
    it = lower_bound(other_segments.begin(), other_segments.end(), t);

    if (it != other_segments.begin()) {
      it--;
      t_final_segment = (it->t_start() < it->t_end() ? it->t_end() : it->t_end() + BETA);
      if (t < t_final_segment) {
        other_length += (t_final_segment < t_final ? t_final_segment - t : t_final - t);
      }
      it++;
    }
    while (it != other_segments.end() && it->t_start() < t_final) {
      t_final_segment = (it->t_start() < it->t_end() ? it->t_end() : it->t_end() + BETA);
      other_length +=
          (t_final_segment < t_final ? t_final_segment - it->t_start() : t_final - it->t_start());
      it++;
    }
    // check if last segment overlaps
    it = other_segments.end();
    it--;
    if (it->t_end() < it->t_start() && t < it->t_end()) {
      other_length += (t_final < it->t_end() ? t_final - t : it->t_end() - t);
    }
  }
  return other_length;
}

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
template <typename orbital_configuration_t>
void ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::compute_intervals(
    double t, double BETA, double& t_up, double& t_down, orbital_configuration_t& segments,
    typename orbital_configuration_t::iterator& s_up,
    typename orbital_configuration_t::iterator& s_down) {
  if (segments.size() == 0) {
    t_up = BETA;
    t_down = BETA;
    s_up = segments.end();
    s_down = segments.end();
  }
  else {
    s_up = lower_bound(segments.begin(), segments.end(), t);

    if (s_up == segments.begin()) {
      s_down = segments.end();
      s_down--;
      if (s_down->t_end() < s_down->t_start())
        t_down = t - s_down->t_end();
      else
        t_down = t + BETA - s_down->t_end();
    }
    else {
      s_down = s_up;
      s_down--;
      if (s_down->t_end() > s_down->t_start())
        t_down = t - s_down->t_end();
      else
        t_down = t - (BETA + s_down->t_end());
    }

    if (s_up == segments.end()) {
      t_up = BETA - t + segments.begin()->t_start();
    }
    else {
      t_up = s_up->t_start() - t;
    }
  }
}

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
template <typename Hybridization_function_t, typename vertex_vertex_matrix_type,
          typename orbital_configuration_t>
double ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::det_rat_up(
    int this_flavor, Hybridization_vertex& new_segment, vertex_vertex_matrix_type& M,
    orbital_configuration_t& segments_old, Hybridization_function_t& F, std::vector<double>& R,
    std::vector<double>& Q_prime, double& det_rat_sign, double& overlap) {
  int inc = 1;
  std::vector<double> Q(1, 0.);

  if (M.size().first > (int)Q.size())
    Q.resize(M.size().first);

  // int* coor = new int[2];
  int coor[2];

  nu nu_obj;

  nu_obj.linind_2_subind(this_flavor, coor);

  // S = F(\tau^n_{e} - \tau^n_{s})
  // R[i] =  F(\tau^n_{e} - \tau^i_{s})
  // Q[i] =  F(\tau^i_{e} - \tau^n_{s})
  double det_rat = interpolate_F(coor, new_segment.t_end() - new_segment.t_start(), F);

  if (M.size().first > 0) {
    typename orbital_configuration_t::iterator it = segments_old.begin();
    for (size_t i = 0; i < segments_old.size(); i++) {
      R[i] = interpolate_F(coor, new_segment.t_end() - it->t_start(), F);
      Q[i] = interpolate_F(coor, it->t_end() - new_segment.t_start(), F);
      it++;
    }

    compute_Q_prime(Q, M, Q_prime);

    int size = M.size().first;
    det_rat += dca::linalg::blas::dot(size, &R[0], inc, &Q_prime[0], inc);
  }

  // take care of sign changes produced by segments which "wind around"
  {
    if (new_segment.t_end() < new_segment.t_start()) {
      det_rat *= -1;
      overlap = -1;
    }
    else
      overlap = 1;

    if (det_rat < 0) {
      det_rat *= -1;
      det_rat_sign = -1;
    }
    else
      det_rat_sign = 1;
  }

  return det_rat;
}

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
template <typename Hybridization_function_t, typename vertex_vertex_matrix_type,
          typename orbital_configuration_t>
void ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::compute_M_up(
    int r, int s, vertex_vertex_matrix_type& M, orbital_configuration_t& /*segments_old*/,
    Hybridization_function_t& /*F*/, std::vector<double>& R, std::vector<double>& Q_prime,
    double S_prime_inv) {
  std::vector<double> R_prime(1, 0.);

  int i_new, j_new;
  int size = M.size().first;

  if (size > 0) {
    if (size > (int)R_prime.size())
      R_prime.resize(M.size().first);

    compute_R_prime(R, M, R_prime);

    compute_M(Q_prime, R_prime, 1. / S_prime_inv, M);
  }

  if (r == 0 && s != 0) {  // need to permute indices of R, L, M
    cycle_column_forward(M);
    // M.cycle_column_forward();
  }

  dca::linalg::matrixop::insertRow(M, r);
  dca::linalg::matrixop::insertCol(M, s);

  // row k+1 and column k
  if (r != 0 || r == s) {  // segments remain in the usual order
    for (int i = 0; i < size; i++) {
      i_new = (i < r ? i : i + 1);
      j_new = (i < s ? i : i + 1);

      M(i_new, s) = Q_prime[i] / S_prime_inv;
      M(r, j_new) = R_prime[i] / S_prime_inv;
    }
  }
  else {  // need to permute indices of R, L, M
    for (int i = 0; i < size; i++) {
      i_new = (i < r ? i : i + 1);
      j_new = (i < s ? i : i + 1);

      M(i_new, s) = Q_prime[i] / S_prime_inv;
      M(r, j_new) =
          R_prime[ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::cycle(
              i, size)] /
          S_prime_inv;
    }
  }

  // fill S_prime
  M(r, s) = 1. / S_prime_inv;
}

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
template <typename vertex_vertex_matrix_type, typename orbital_configuration_t>
double ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::det_rat_down(
    int r, int s, vertex_vertex_matrix_type& M, orbital_configuration_t& segments_old,
    double& det_rat_sign) {
  double det_rat = M(r, s);

  // take care of sign changes produced by segments which "wind around"
  if (r == int(segments_old.size()) - 1) {
    typename orbital_configuration_t::iterator it = segments_old.end();
    it--;
    if (it->t_end() < it->t_start())
      det_rat *= -1;
  }

  if (det_rat < 0) {
    det_rat_sign = -1;
    det_rat *= -1;
  }
  else {
    det_rat_sign = 1;
  }

  return det_rat;
}

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
template <typename vertex_vertex_matrix_type>
void ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::compute_M_down(
    int r, int s, vertex_vertex_matrix_type& M) {
  int inc = 1;

  std::vector<double> Q_prime(1, 0.);
  std::vector<double> R_prime(1, 0.);

  if (M.size().first > 1) {
    if (M.size().first > (int)R_prime.size()) {
      Q_prime.resize(M.size().first);
      R_prime.resize(M.size().first);
    }

    int incy = M.leadingDimension();
    int size = M.size().first;

    dca::linalg::blas::copy(size, &M(0, s), inc, &Q_prime[0], inc);
    dca::linalg::blas::copy(size, &M(r, 0), incy, &R_prime[0], inc);

    compute_M(Q_prime, R_prime, -1. / Q_prime[r], M);

    dca::linalg::matrixop::removeRowAndCol(M, r, s);

    if (r == 0 && s != 0) {  // need to permute indices of R, L, M
      // M.cycle_column_backward();
      cycle_column_backward(M);
    }
  }
  else
    M.resize(0);
}

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
template <typename Hybridization_function_t, typename vertex_vertex_matrix_type,
          typename orbital_configuration_t>
double ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::det_rat_shift_end(
    int this_flavor, Hybridization_vertex& new_segment, int k, vertex_vertex_matrix_type& M,
    orbital_configuration_t& segments_old, Hybridization_function_t& F, std::vector<double>& R,
    std::vector<double>& Q_prime, double& det_rat_sign, double& overlap) {
  int inc = 1;
  int coor[2];

  nu nu_obj;
  nu_obj.linind_2_subind(this_flavor, coor);

  assert(M.size().first == M.size().second);
  int size = M.size().first;

  dca::linalg::blas::copy(size, &M(0, k), inc, &Q_prime[0], inc);

  typename orbital_configuration_t::iterator it0;
  it0 = segments_old.begin();

  for (int i = 0; i < M.size().first; i++) {
    R[i] = interpolate_F(coor, new_segment.t_end() - (it0 + i)->t_start(), F);
  }

  double det_rat = dca::linalg::blas::dot(size, &R[0], inc, &Q_prime[0], inc);

  for (int i = 0; i < M.size().first; i++) {
    R[i] = interpolate_F(coor, new_segment.t_end() - (it0 + i)->t_start(), F) -
           interpolate_F(coor, (it0 + k)->t_end() - (it0 + i)->t_start(), F);
  }

  {
    overlap = 1;
    //  take care of sign changes produced by segments which "wind around"
    if ((new_segment.t_end() - new_segment.t_start()) * ((it0 + k)->t_end() - (it0 + k)->t_start()) <
        0) {
      det_rat *= -1;
      overlap = -1;
    }

    if (det_rat < 0) {
      det_rat_sign = -1;
      det_rat *= -1;
    }
    else
      det_rat_sign = 1;
  }

  return det_rat;
}

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
template <typename Hybridization_function_t, typename vertex_vertex_matrix_type,
          typename orbital_configuration_t>
double ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::det_rat_shift_start(
    int this_flavor, Hybridization_vertex& new_segment, int k, vertex_vertex_matrix_type& M,
    orbital_configuration_t& segments_old, Hybridization_function_t& F,
    std::vector<double>& R_prime, std::vector<double>& Q, double& det_rat_sign, double& overlap) {
  int inc = 1;
  int coor[2];

  nu nu_obj;
  nu_obj.linind_2_subind(this_flavor, coor);

  assert(M.size().first == M.size().second);
  int size = M.size().first;
  int incy = M.leadingDimension();

  dca::linalg::blas::copy(size, &M(k, 0), incy, &R_prime[0], inc);

  typename orbital_configuration_t::iterator it0;
  it0 = segments_old.begin();

  for (int i = 0; i < M.size().first; i++) {
    Q[i] = interpolate_F(coor, (it0 + i)->t_end() - new_segment.t_start(), F);
  }

  double det_rat = dca::linalg::blas::dot(size, &R_prime[0], inc, &Q[0], inc);

  for (int i = 0; i < M.size().first; i++) {
    Q[i] = interpolate_F(coor, (it0 + i)->t_end() - new_segment.t_start(), F) -
           interpolate_F(coor, (it0 + i)->t_end() - (it0 + k)->t_start(), F);
  }

  {
    overlap = 1;
    // take care of sign changes produced by segments which "wind around"

    if ((new_segment.t_end() - new_segment.t_start()) * ((it0 + k)->t_end() - (it0 + k)->t_start()) <
        0) {
      det_rat *= -1;
      overlap = -1;
    }
    // std::cout<<"\n\t"<<det_rat<<" "<<" "<<"\n";

    if (det_rat < 0) {
      det_rat_sign = -1;
      det_rat *= -1;
    }
    else
      det_rat_sign = 1;
  }

  return det_rat;
}

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
template <typename vertex_vertex_matrix_type>
void ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::compute_M_shift_end(
    int /*k*/, vertex_vertex_matrix_type& M, std::vector<double>& R, std::vector<double>& Q_prime,
    double det_rat) {
  // int inc = 1;

  std::vector<double> R_prime(1, 0.);

  double S_prime = 1. / det_rat;

  if (M.size().first > (int)R_prime.size())
    R_prime.resize(M.size().first);

  compute_R_prime(R, M, R_prime);

  compute_M(Q_prime, R_prime, S_prime, M);
}

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
template <typename vertex_vertex_matrix_type>
void ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::compute_M_shift_start(
    int /*k*/, vertex_vertex_matrix_type& M, std::vector<double>& R_prime, std::vector<double>& Q,
    double det_rat) {
  // int inc = 1;

  std::vector<double> Q_prime(1, 0.);

  double S_prime = 1. / det_rat;

  if (M.size().first > (int)Q_prime.size())
    Q_prime.resize(M.size().first);

  compute_Q_prime(Q, M, Q_prime);

  compute_M(Q_prime, R_prime, S_prime, M);
}

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
template <typename vertex_vertex_matrix_type>
void ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::compute_Q_prime(
    std::vector<double>& Q, vertex_vertex_matrix_type& M, std::vector<double>& Q_prime) {
  dca::linalg::blas::gemv("N", M.size().first, M.size().second, -1., &M(0, 0), M.leadingDimension(),
                          &Q[0], 1, 0., &Q_prime[0], 1);
}

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
template <typename vertex_vertex_matrix_type>
void ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::compute_R_prime(
    std::vector<double>& R, vertex_vertex_matrix_type& M, std::vector<double>& R_prime) {
  dca::linalg::blas::gemv("T", M.size().first, M.size().second, -1., &M(0, 0), M.leadingDimension(),
                          &R[0], 1, 0., &R_prime[0], 1);
}

template <typename parameters_t, typename MOMS_t, typename configuration_t, typename rng_t>
template <typename vertex_vertex_matrix_type>
void ss_hybridization_walker_routines<parameters_t, MOMS_t, configuration_t, rng_t>::compute_M(
    std::vector<double>& Q_prime, std::vector<double>& R_prime, double S_prime,
    vertex_vertex_matrix_type& M) {
  dca::linalg::blas::gemm("N", "N", M.size().first, M.size().second, 1, S_prime, &Q_prime[0],
                          M.leadingDimension(), &R_prime[0], 1, 1., &M(0, 0), M.leadingDimension());
  // dca::linalg::blas::ger(M.size().first, M.size().second, S_prime,
  // &Q_prime[0], 1, &R_prime[0], 1, &M(0, 0), M.leadingDimension());
}

}  // cthyb
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SS_CT_HYB_WALKER_TOOLS_SS_HYBRIDIZATION_WALKER_ROUTINES_HPP
