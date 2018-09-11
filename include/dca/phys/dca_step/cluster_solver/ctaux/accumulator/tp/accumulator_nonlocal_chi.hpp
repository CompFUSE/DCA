// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class computes the nonlocal \f$\chi(k_1,k_2,q)\f$.

/*
 *  Definition of two-particle functions:
 *
 *  \section ph particle-hole channel:
 *
 *  The magnetic channel,
 *
 *  \f{eqnarray*}{
 *    G^{II}_{ph} &=& \frac{1}{4} \sum_{\sigma_1,\sigma_2 = \pm1} (\sigma_1 \: \sigma_2) \langle
 *    T_\tau\{c^\dagger_{k_1+q,\sigma_1} c_{k_1,\sigma_1} c^\dagger_{k_2,\sigma_2}
 *    c_{k_2+q,\sigma_2} \}\rangle
 *  \f}
 *
 *  The charge channel,
 *
 *  \f{eqnarray*}{
 *    G^{II}_{ph} &=& \frac{1}{4} \sum_{\sigma_1,\sigma_2 = \pm1} \langle
 *    T_\tau\{c^\dagger_{k_1+q,\sigma_1} c_{k_1,\sigma_1} c^\dagger_{k_2,\sigma_2}
 *    c_{k_2+q,\sigma_2} \}\rangle
 *  \f}
 *
 *  The transverse channel,
 *
 *  \f{eqnarray*}{
 *    G^{II}_{ph} &=& \frac{1}{2} \sum_{\sigma = \pm1} \langle T_\tau\{c^\dagger_{k_1+q,\sigma}
 *    c_{k_1,-\sigma} c^\dagger_{k_2,-\sigma}  c_{k_2+q,\sigma} \}\rangle
 *  \f}
 *
 *  \section pp particle-hole channel:
 *
 *  The transverse (or superconducting) channel,
 *
 *  \f{eqnarray*}{
 *    G^{II}_{pp} &=& \frac{1}{2} \sum_{\sigma= \pm1} \langle T_\tau\{ c^\dagger_{q-k_1,\sigma}
 *    c^\dagger_{k_1,-\sigma} c_{k_2,-\sigma}  c_{q-k_2,\sigma} \}
 *  \f}
 */

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_ACCUMULATOR_TP_ACCUMULATOR_NONLOCAL_CHI_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_ACCUMULATOR_TP_ACCUMULATOR_NONLOCAL_CHI_HPP

#include <cassert>
#include <cmath>
#include <complex>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"
#include "dca/phys/four_point_type.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

template <class parameters_type, class MOMS_type>
class accumulator_nonlocal_chi {
public:
  using TpGreensFunction = typename MOMS_type::TpGreensFunction;

  using w_VERTEX = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT>>;
  using w_VERTEX_EXTENDED = func::dmn_0<domains::vertex_frequency_domain<domains::EXTENDED>>;
  using w_VERTEX_EXTENDED_POS =
      func::dmn_0<domains::vertex_frequency_domain<domains::EXTENDED_POSITIVE>>;

  using b = func::dmn_0<domains::electron_band_domain>;

  using CDA = ClusterDomainAliases<parameters_type::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using KClusterDmn = typename CDA::KClusterDmn;

  // This needs shifted to new aliases
  typedef RClusterDmn r_dmn_t;
  typedef KClusterDmn k_dmn_t;
  typedef typename r_dmn_t::parameter_type r_cluster_type;
  typedef typename k_dmn_t::parameter_type k_cluster_type;

  typedef typename parameters_type::profiler_type profiler_t;
  typedef typename parameters_type::concurrency_type concurrency_type;

  typedef typename parameters_type::MC_measurement_scalar_type scalar_type;

  typedef typename parameters_type::G4_w1_dmn_t w1_dmn_t;
  typedef typename parameters_type::G4_w2_dmn_t w2_dmn_t;

  typedef func::dmn_variadic<b, b, r_dmn_t, r_dmn_t, w1_dmn_t, w2_dmn_t> b_b_r_r_w_w_dmn_t;
  typedef func::dmn_variadic<b, b, k_dmn_t, k_dmn_t, w1_dmn_t, w2_dmn_t> b_b_k_k_w_w_dmn_t;

public:
  accumulator_nonlocal_chi(parameters_type& parameters_ref, MOMS_type& MOMS_ref, int id,
                           TpGreensFunction& G4_ref);

  void initialize();

  void finalize();

  template <class nonlocal_G_t>
  void execute(scalar_type current_sign, nonlocal_G_t& nonlocal_G_obj);

private:
  void F(
      int n1, int m1, int k1, int k2, int w1, int w2,
      func::function<std::complex<scalar_type>,
                     func::dmn_variadic<b, b, k_dmn_t, k_dmn_t, w_VERTEX_EXTENDED, w_VERTEX_EXTENDED>>& G2,
      std::complex<scalar_type>& G2_result);

  void F(int n1, int m1, int k1, int k2, int w1, int w2,
         func::function<
             std::complex<scalar_type>,
             func::dmn_variadic<b, b, k_dmn_t, k_dmn_t, w_VERTEX_EXTENDED_POS, w_VERTEX_EXTENDED>>& G2,
         std::complex<scalar_type>& G2_result);

  void F(
      int n1, int m1, int k1, int k2, int w1, int w2,
      func::function<std::complex<scalar_type>,
                     func::dmn_variadic<b, b, k_dmn_t, k_dmn_t, w_VERTEX_EXTENDED, w_VERTEX_EXTENDED>>& G2_dn,
      std::complex<scalar_type>& G2_dn_result,
      func::function<std::complex<scalar_type>,
                     func::dmn_variadic<b, b, k_dmn_t, k_dmn_t, w_VERTEX_EXTENDED, w_VERTEX_EXTENDED>>& G2_up,
      std::complex<scalar_type>& G2_up_result);

  void F(int n1, int m1, int k1, int k2, int w1, int w2,
         func::function<
             std::complex<scalar_type>,
             func::dmn_variadic<b, b, k_dmn_t, k_dmn_t, w_VERTEX_EXTENDED_POS, w_VERTEX_EXTENDED>>& G2_dn,
         std::complex<scalar_type>& G2_dn_result,
         func::function<
             std::complex<scalar_type>,
             func::dmn_variadic<b, b, k_dmn_t, k_dmn_t, w_VERTEX_EXTENDED_POS, w_VERTEX_EXTENDED>>& G2_up,
         std::complex<scalar_type>& G2_up_result);

  void accumulate_particle_hole_transverse(
      func::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& G2_k_k_w_w_e_DN,
      func::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& G2_k_k_w_w_e_UP,
      scalar_type sign);

  void accumulate_particle_hole_magnetic(
      func::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& G2_k_k_w_w_e_DN,
      func::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& G2_k_k_w_w_e_UP,
      scalar_type sign);

  void accumulate_particle_hole_charge(
      func::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& G2_k_k_w_w_e_DN,
      func::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& G2_k_k_w_w_e_UP,
      scalar_type sign);

  void accumulate_particle_particle_superconducting(
      func::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& G2_k_k_w_w_e_DN,
      func::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& G2_k_k_w_w_e_UP,
      scalar_type sign);

private:
  parameters_type& parameters;
  MOMS_type& MOMS;
  concurrency_type& concurrency;

  int thread_id;

  TpGreensFunction& G4;

  int w_VERTEX_EXTENDED_POS_dmn_size;

  func::dmn_variadic<b, b, k_dmn_t, k_dmn_t, w_VERTEX_EXTENDED, w_VERTEX_EXTENDED> b_b_k_k_w_full_w_full_dmn;
  func::dmn_variadic<b, b, k_dmn_t, k_dmn_t, w_VERTEX_EXTENDED_POS, w_VERTEX_EXTENDED> b_b_k_k_w_pos_w_full_dmn;

  func::function<int, k_dmn_t> min_k_dmn_t;
  func::function<int, k_dmn_t> q_plus_;
  func::function<int, k_dmn_t> q_min_;

  func::function<int, w_VERTEX> min_w_vertex;
  func::function<int, w_VERTEX_EXTENDED> min_w_vertex_ext;

  func::function<int, w_VERTEX> w_vertex_2_w_vertex_ext;

  func::function<int, w_VERTEX_EXTENDED> w_vertex_ext_2_w_vertex_ext_pos;
};

template <class parameters_type, class MOMS_type>
accumulator_nonlocal_chi<parameters_type, MOMS_type>::accumulator_nonlocal_chi(
    parameters_type& parameters_ref, MOMS_type& MOMS_ref, int id, TpGreensFunction& G4_ref)
    : parameters(parameters_ref),
      MOMS(MOMS_ref),
      concurrency(parameters.get_concurrency()),

      thread_id(id),

      G4(G4_ref),

      w_VERTEX_EXTENDED_POS_dmn_size(w_VERTEX_EXTENDED_POS::dmn_size()),

      b_b_k_k_w_full_w_full_dmn(),
      b_b_k_k_w_pos_w_full_dmn(),

      min_k_dmn_t("min_k_dmn_t"),
      q_plus_("q_plus_"),
      q_min_("q_min_"),

      min_w_vertex(" min_w_vertex"),
      min_w_vertex_ext("min_w_vertex_ext"),

      w_vertex_2_w_vertex_ext("w_vertex_2_w_vertex_ext"),

      w_vertex_ext_2_w_vertex_ext_pos("w_vertex_ext_2_w_vertex_ext_pos") {
  int q_channel = parameters.get_four_point_momentum_transfer_index();
  // int k0_index  = k_cluster_type::get_k_0_index();
  int k0_index = k_cluster_type::origin_index();

  for (int l = 0; l < k_cluster_type::get_size(); l++) {
    min_k_dmn_t(l) = k_cluster_type::subtract(l, k0_index);

    q_plus_(l) = k_cluster_type::add(l, q_channel);
    q_min_(l) = k_cluster_type::subtract(l, q_channel);
  }

  {
    for (int l = 0; l < w_VERTEX::dmn_size(); l++)
      min_w_vertex(l) = w_VERTEX::dmn_size() - 1 - l;

    for (int l = 0; l < w_VERTEX_EXTENDED::dmn_size(); l++)
      min_w_vertex_ext(l) = w_VERTEX_EXTENDED::dmn_size() - 1 - l;
  }

  {
    for (int i = 0; i < w_VERTEX::dmn_size(); i++)
      for (int j = 0; j < w_VERTEX_EXTENDED::dmn_size(); j++)
        if (std::fabs(w_VERTEX::get_elements()[i] - w_VERTEX_EXTENDED::get_elements()[j]) < 1.e-6)
          w_vertex_2_w_vertex_ext(i) = j;
  }

  {
    for (int l = 0; l < w_VERTEX_EXTENDED::dmn_size(); l++) {
      if (l < w_VERTEX_EXTENDED_POS::dmn_size())
        w_vertex_ext_2_w_vertex_ext_pos(l) = w_VERTEX_EXTENDED_POS::dmn_size() - 1 - l;
      else
        w_vertex_ext_2_w_vertex_ext_pos(l) = l - w_VERTEX_EXTENDED_POS::dmn_size();

      // cout << l << "\t" << w_vertex_ext_2_w_vertex_ext_pos(l) << "\n";
    }
  }
}

template <class parameters_type, class MOMS_type>
void accumulator_nonlocal_chi<parameters_type, MOMS_type>::initialize() {
  MOMS.get_G4() = 0.;
}

template <class parameters_type, class MOMS_type>
void accumulator_nonlocal_chi<parameters_type, MOMS_type>::finalize() {
  //       for(int i=0; i<G4.size(); i++)
  //      MOMS.G4_k_k_w_w(i) = G4(i);
}

template <class parameters_type, class MOMS_type>
template <class nonlocal_G_t>
void accumulator_nonlocal_chi<parameters_type, MOMS_type>::execute(scalar_type current_sign,
                                                                   nonlocal_G_t& nonlocal_G_obj) {
  profiler_t profiler("compute nonlocal-chi", "CT-AUX accumulator", __LINE__, thread_id);

  switch (parameters.get_four_point_type()) {
    case PARTICLE_HOLE_TRANSVERSE:
      accumulate_particle_hole_transverse(nonlocal_G_obj.get_G_k_k_w_w_e_DN(),
                                          nonlocal_G_obj.get_G_k_k_w_w_e_UP(), current_sign);
      break;

    case PARTICLE_HOLE_MAGNETIC:
      accumulate_particle_hole_magnetic(nonlocal_G_obj.get_G_k_k_w_w_e_DN(),
                                        nonlocal_G_obj.get_G_k_k_w_w_e_UP(), current_sign);
      break;

    case PARTICLE_HOLE_CHARGE:
      accumulate_particle_hole_charge(nonlocal_G_obj.get_G_k_k_w_w_e_DN(),
                                      nonlocal_G_obj.get_G_k_k_w_w_e_UP(), current_sign);
      break;

    case PARTICLE_PARTICLE_UP_DOWN:
      accumulate_particle_particle_superconducting(
          nonlocal_G_obj.get_G_k_k_w_w_e_DN(), nonlocal_G_obj.get_G_k_k_w_w_e_UP(), current_sign);
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

template <class parameters_type, class MOMS_type>
inline void accumulator_nonlocal_chi<parameters_type, MOMS_type>::F(
    int n1, int m1, int k1, int k2, int w1, int w2,
    func::function<std::complex<scalar_type>,
                   func::dmn_variadic<b, b, k_dmn_t, k_dmn_t, w_VERTEX_EXTENDED, w_VERTEX_EXTENDED>>& G2,
    std::complex<scalar_type>& G2_result) {
  G2_result = G2(n1, m1, k1, k2, w1, w2);
}

template <class parameters_type, class MOMS_type>
inline void accumulator_nonlocal_chi<parameters_type, MOMS_type>::F(
    int n1, int m1, int k1, int k2, int w1, int w2,
    func::function<std::complex<scalar_type>,
                   func::dmn_variadic<b, b, k_dmn_t, k_dmn_t, w_VERTEX_EXTENDED_POS, w_VERTEX_EXTENDED>>& G2,
    std::complex<scalar_type>& G2_result) {
  if (w1 < w_VERTEX_EXTENDED_POS_dmn_size) {
    G2_result = conj(G2(n1, m1, min_k_dmn_t(k1), min_k_dmn_t(k2),
                        w_vertex_ext_2_w_vertex_ext_pos(w1), min_w_vertex_ext(w2)));
  }
  else {
    G2_result = G2(n1, m1, k1, k2, w_vertex_ext_2_w_vertex_ext_pos(w1), w2);
  }
}

template <class parameters_type, class MOMS_type>
inline void accumulator_nonlocal_chi<parameters_type, MOMS_type>::F(
    int n1, int m1, int k1, int k2, int w1, int w2,
    func::function<std::complex<scalar_type>,
                   func::dmn_variadic<b, b, k_dmn_t, k_dmn_t, w_VERTEX_EXTENDED, w_VERTEX_EXTENDED>>& G2_dn,
    std::complex<scalar_type>& G2_dn_result,
    func::function<std::complex<scalar_type>,
                   func::dmn_variadic<b, b, k_dmn_t, k_dmn_t, w_VERTEX_EXTENDED, w_VERTEX_EXTENDED>>& G2_up,
    std::complex<scalar_type>& G2_up_result) {
  int lin_ind = b_b_k_k_w_full_w_full_dmn(n1, m1, k1, k2, w1, w2);

  G2_dn_result = G2_dn(lin_ind);
  G2_up_result = G2_up(lin_ind);
}

template <class parameters_type, class MOMS_type>
inline void accumulator_nonlocal_chi<parameters_type, MOMS_type>::F(
    int n1, int m1, int k1, int k2, int w1, int w2,
    func::function<std::complex<scalar_type>,
                   func::dmn_variadic<b, b, k_dmn_t, k_dmn_t, w_VERTEX_EXTENDED_POS, w_VERTEX_EXTENDED>>& G2_dn,
    std::complex<scalar_type>& G2_dn_result,
    func::function<std::complex<scalar_type>,
                   func::dmn_variadic<b, b, k_dmn_t, k_dmn_t, w_VERTEX_EXTENDED_POS, w_VERTEX_EXTENDED>>& G2_up,
    std::complex<scalar_type>& G2_up_result) {
  if (w1 < w_VERTEX_EXTENDED_POS_dmn_size) {
    int lin_ind = b_b_k_k_w_pos_w_full_dmn(n1, m1, min_k_dmn_t(k1), min_k_dmn_t(k2),
                                           w_vertex_ext_2_w_vertex_ext_pos(w1), min_w_vertex_ext(w2));

    G2_dn_result = conj(G2_dn(lin_ind));
    G2_up_result = conj(G2_up(lin_ind));
  }
  else {
    int lin_ind = b_b_k_k_w_pos_w_full_dmn(n1, m1, k1, k2, w_vertex_ext_2_w_vertex_ext_pos(w1), w2);

    G2_dn_result = G2_dn(lin_ind);
    G2_up_result = G2_up(lin_ind);
  }
}

template <class parameters_type, class MOMS_type>
void accumulator_nonlocal_chi<parameters_type, MOMS_type>::accumulate_particle_hole_magnetic(
    func::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& G2_k_k_w_w_e_DN,
    func::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& G2_k_k_w_w_e_UP, scalar_type sign) {
  // n1 ------------------------ m1
  //        |           |
  //        |           |
  //        |           |
  // n2 ------------------------ m2

  std::complex<scalar_type> G2_DN_n1_m2_k1_k2_w1_w2, G2_UP_n1_m2_k1_k2_w1_w2,
      G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1, G2_UP_n2_m1_k2_plus_q_k1_plus_q_w2_w1,
      G2_DN_n1_m1_k1_k1_plus_q_w1_w1, G2_UP_n1_m1_k1_k1_plus_q_w1_w1,
      G2_DN_n2_m2_k2_plus_q_k2_w2_w2, G2_UP_n2_m2_k2_plus_q_k2_w2_w2, G4_val;

  int k2_plus_q, k1_plus_q;

  scalar_type sign_div_2 = scalar_type(sign) / 2.;

  int w_nu = parameters.get_four_point_frequency_transfer();

  for (int w2 = 0; w2 < w_VERTEX::dmn_size(); w2++) {
    int w2_ext = w_vertex_2_w_vertex_ext(w2);
    int w2_ext_plus_w_nu = w2_ext + w_nu;
    assert(std::fabs(w_VERTEX::get_elements()[w2] - w_VERTEX_EXTENDED::get_elements()[w2_ext]) <
           1.e-6);

    for (int w1 = 0; w1 < w_VERTEX::dmn_size(); w1++) {
      int w1_ext = w_vertex_2_w_vertex_ext(w1);
      int w1_ext_plus_w_nu = w1_ext + w_nu;

      for (int k2 = 0; k2 < k_dmn_t::dmn_size(); k2++) {
        //          int k2_plus_q = k_cluster_type::add(k2,q_channel);
        k2_plus_q = q_plus_(k2);

        for (int k1 = 0; k1 < k_dmn_t::dmn_size(); k1++) {
          //          int k1_plus_q = k_cluster_type::add(k1,q_channel);
          k1_plus_q = q_plus_(k1);

          for (int m1 = 0; m1 < b::dmn_size(); m1++) {
            for (int m2 = 0; m2 < b::dmn_size(); m2++) {
              for (int n1 = 0; n1 < b::dmn_size(); n1++) {
                for (int n2 = 0; n2 < b::dmn_size(); n2++) {
                  F(n1, m2, k1, k2, w1_ext, w2_ext, G2_k_k_w_w_e_DN, G2_DN_n1_m2_k1_k2_w1_w2,
                    G2_k_k_w_w_e_UP, G2_UP_n1_m2_k1_k2_w1_w2);

                  F(n2, m1, k2_plus_q, k1_plus_q, w2_ext_plus_w_nu, w1_ext_plus_w_nu,
                    G2_k_k_w_w_e_DN, G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1, G2_k_k_w_w_e_UP,
                    G2_UP_n2_m1_k2_plus_q_k1_plus_q_w2_w1);

                  F(n1, m1, k1, k1_plus_q, w1_ext, w1_ext_plus_w_nu, G2_k_k_w_w_e_DN,
                    G2_DN_n1_m1_k1_k1_plus_q_w1_w1, G2_k_k_w_w_e_UP, G2_UP_n1_m1_k1_k1_plus_q_w1_w1);

                  F(n2, m2, k2_plus_q, k2, w2_ext_plus_w_nu, w2_ext, G2_k_k_w_w_e_DN,
                    G2_DN_n2_m2_k2_plus_q_k2_w2_w2, G2_k_k_w_w_e_UP, G2_UP_n2_m2_k2_plus_q_k2_w2_w2);

                  G4_val = -(G2_DN_n1_m2_k1_k2_w1_w2 * G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1 +
                             G2_UP_n1_m2_k1_k2_w1_w2 * G2_UP_n2_m1_k2_plus_q_k1_plus_q_w2_w1)

                           +
                           (G2_UP_n1_m1_k1_k1_plus_q_w1_w1 - G2_DN_n1_m1_k1_k1_plus_q_w1_w1) *
                               (G2_UP_n2_m2_k2_plus_q_k2_w2_w2 - G2_DN_n2_m2_k2_plus_q_k2_w2_w2);

                  /*
                    G4 = - (G2_k_k_w_w_e_DN(n1, m2, k1, k2, w1, w2) * G2_k_k_w_w_e_DN(n2, m1,
                    k2_plus_q, k1_plus_q, w2+w_channel, w1+w_channel)
                    + G2_k_k_w_w_e_UP(n1, m2, k1, k2, w1, w2) * G2_k_k_w_w_e_UP(n2, m1, k2_plus_q,
                    k1_plus_q, w2+w_channel, w1+w_channel))

                    + (G2_k_k_w_w_e_UP(n1, m1, k1, k1_plus_q, w1, w1+w_channel) -
                    G2_k_k_w_w_e_DN(n1, m1, k1, k1_plus_q, w1, w1+w_channel))
                    * (G2_k_k_w_w_e_UP(n2, m2, k2_plus_q, k2, w2+w_channel, w2) -
                    G2_k_k_w_w_e_DN(n2, m2, k2_plus_q, k2, w2+w_channel, w2));
                    */

                  G4(n1, n2, m1, m2, k1, k2, 0, w1, w2, 0) +=
                      std::complex<double>(sign_div_2 * G4_val);
                  // MOMS.G4_k_k_w_w(n1, n2, m1, m2, k1, k2, w1, w2) +=
                  // std::complex<double>(sign_div_2 * G4);
                }
              }
            }
          }
        }
      }
    }
  }
}

template <class parameters_type, class MOMS_type>
void accumulator_nonlocal_chi<parameters_type, MOMS_type>::accumulate_particle_hole_transverse(
    func::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& G2_k_k_w_w_e_DN,
    func::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& G2_k_k_w_w_e_UP, scalar_type sign) {
  // n1 ------------------------ m1
  //        |           |
  //        |           |
  //        |           |
  // n2 ------------------------ m2

  std::complex<scalar_type> G2_DN_n1_m2_k1_k2_w1_w2, G2_UP_n1_m2_k1_k2_w1_w2,
      G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1, G2_UP_n2_m1_k2_plus_q_k1_plus_q_w2_w1, G4_val;

  int k2_plus_q, k1_plus_q;

  scalar_type sign_div_2 = scalar_type(sign) / 2.;

  int w_nu = parameters.get_four_point_frequency_transfer();

  for (int w2 = 0; w2 < w_VERTEX::dmn_size(); w2++) {
    int w2_ext = w_vertex_2_w_vertex_ext(w2);
    int w2_ext_plus_w_nu = w2_ext + w_nu;
    assert(std::fabs(w_VERTEX::get_elements()[w2] - w_VERTEX_EXTENDED::get_elements()[w2_ext]) <
           1.e-6);

    for (int w1 = 0; w1 < w_VERTEX::dmn_size(); w1++) {
      int w1_ext = w_vertex_2_w_vertex_ext(w1);
      int w1_ext_plus_w_nu = w1_ext + w_nu;

      for (int k2 = 0; k2 < k_dmn_t::dmn_size(); k2++) {
        k2_plus_q = q_plus_(k2);

        for (int k1 = 0; k1 < k_dmn_t::dmn_size(); k1++) {
          k1_plus_q = q_plus_(k1);

          for (int m1 = 0; m1 < b::dmn_size(); m1++) {
            for (int m2 = 0; m2 < b::dmn_size(); m2++) {
              for (int n1 = 0; n1 < b::dmn_size(); n1++) {
                for (int n2 = 0; n2 < b::dmn_size(); n2++) {
                  //                    F(n1, m2, k1, k2, w1, w2,
                  //                      G2_k_k_w_w_e_DN, G2_DN_n1_m2_k1_k2_w1_w2,
                  //                      G2_k_k_w_w_e_UP, G2_UP_n1_m2_k1_k2_w1_w2);

                  //                    F(n2, m1, k2_plus_q, k1_plus_q, w2+w_nu, w1+w_nu,
                  //                      G2_k_k_w_w_e_DN, G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1,
                  //                      G2_k_k_w_w_e_UP, G2_UP_n2_m1_k2_plus_q_k1_plus_q_w2_w1);

                  F(n1, m2, k1, k2, w1_ext, w2_ext, G2_k_k_w_w_e_DN, G2_DN_n1_m2_k1_k2_w1_w2,
                    G2_k_k_w_w_e_UP, G2_UP_n1_m2_k1_k2_w1_w2);

                  F(n2, m1, k2_plus_q, k1_plus_q, w2_ext_plus_w_nu, w1_ext_plus_w_nu,
                    G2_k_k_w_w_e_DN, G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1, G2_k_k_w_w_e_UP,
                    G2_UP_n2_m1_k2_plus_q_k1_plus_q_w2_w1);

                  G4_val = -(G2_DN_n1_m2_k1_k2_w1_w2 * G2_UP_n2_m1_k2_plus_q_k1_plus_q_w2_w1 +
                             G2_UP_n1_m2_k1_k2_w1_w2 * G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1);

                  G4(n1, n2, m1, m2, k1, k2, 0, w1, w2, 0) +=
                      std::complex<double>(sign_div_2 * G4_val);
                }
              }
            }
          }
        }
      }
    }
  }
}

/*
  template<class parameters_type, class MOMS_type>
  void accumulator_nonlocal_chi<parameters_type,
  MOMS_type>::accumulate_particle_hole_magnetic_fast(func::function<std::complex<scalar_type>,
  func::dmn_variadic<b,b,k_dmn_t,k_dmn_t,w_VERTEX,w_VERTEX> >& G2_k_k_w_w_e_DN,
  func::function<std::complex<scalar_type>,
  func::dmn_variadic<b,b,k_dmn_t,k_dmn_t,w_VERTEX,w_VERTEX> >&
  G2_k_k_w_w_e_UP,
  scalar_type sign)
  {
  std::complex<scalar_type> *G2_DN_n1_m2_k1_k2_w1_w2, *G2_UP_n1_m2_k1_k2_w1_w2,
  *G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1, *G2_UP_n2_m1_k2_plus_q_k1_plus_q_w2_w1,
  *G2_DN_n1_m1_k1_k1_plus_q_w1_w1, *G2_UP_n1_m1_k1_k1_plus_q_w1_w1,
  *G2_DN_n2_m2_k2_plus_q_k2_w2_w2, *G2_UP_n2_m2_k2_plus_q_k2_w2_w2;

  std::complex<double>* G4;

  int k2_plus_q, k1_plus_q;

  scalar_type sign_div_2 =scalar_type(sign)/2.;

  for(int w2=0; w2<w_VERTEX::dmn_size(); w2++){
  for(int w1=0; w1<w_VERTEX::dmn_size(); w1++){

  for(int k2=0; k2<k_dmn_t::dmn_size(); k2++){

  //        int k2_plus_q = k_cluster_type::add(k2,q_channel);
  k2_plus_q = q_plus_(k2);

  for(int k1=0; k1<k_dmn_t::dmn_size(); k1++){

  //          int k1_plus_q = k_cluster_type::add(k1,q_channel);
  k1_plus_q = q_plus_(k1);

  G2_DN_n1_m2_k1_k2_w1_w2 = &G2_k_k_w_w_e_DN(0, 0, k1, k2, w1, w2);
  G2_UP_n1_m2_k1_k2_w1_w2 = &G2_k_k_w_w_e_UP(0, 0, k1, k2, w1, w2);

  G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1 = &G2_k_k_w_w_e_DN(0, 0, k2_plus_q, k1_plus_q, w2, w1);
  G2_UP_n2_m1_k2_plus_q_k1_plus_q_w2_w1 = &G2_k_k_w_w_e_UP(0, 0, k2_plus_q, k1_plus_q, w2, w1);

  G2_DN_n1_m1_k1_k1_plus_q_w1_w1 = &G2_k_k_w_w_e_DN(0, 0, k1, k1_plus_q, w1, w1);
  G2_UP_n1_m1_k1_k1_plus_q_w1_w1 = &G2_k_k_w_w_e_UP(0, 0, k1, k1_plus_q, w1, w1);

  G2_DN_n2_m2_k2_plus_q_k2_w2_w2 = &G2_k_k_w_w_e_DN(0, 0, k2_plus_q, k2, w2, w2);
  G2_UP_n2_m2_k2_plus_q_k2_w2_w2 = &G2_k_k_w_w_e_UP(0, 0, k2_plus_q, k2, w2, w2);

  G4 = &MOMS.G4_k_k_w_w(0,0,0,0, k1, k2, w1, w2);

  accumulator_nonlocal_chi_atomic<model, PARTICLE_HOLE_MAGNETIC>::execute(G2_DN_n1_m2_k1_k2_w1_w2
  , G2_UP_n1_m2_k1_k2_w1_w2,
  G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1, G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1,
  G2_DN_n1_m1_k1_k1_plus_q_w1_w1       , G2_UP_n1_m1_k1_k1_plus_q_w1_w1,
  G2_DN_n2_m2_k2_plus_q_k2_w2_w2       , G2_UP_n2_m2_k2_plus_q_k2_w2_w2,
  G4, sign_div_2);
  }
  }
  }
  }
  }
*/

template <class parameters_type, class MOMS_type>
void accumulator_nonlocal_chi<parameters_type, MOMS_type>::accumulate_particle_hole_charge(
    func::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& G2_k_k_w_w_e_DN,
    func::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& G2_k_k_w_w_e_UP, scalar_type sign) {
  // n1 ------------------------ m1
  //        |           |
  //        |           |
  //        |           |
  // n2 ------------------------ m2

  // int q_channel = parameters.get_q_channel();

  std::complex<scalar_type> G2_DN_n1_m2_k1_k2_w1_w2, G2_UP_n1_m2_k1_k2_w1_w2,
      G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1, G2_UP_n2_m1_k2_plus_q_k1_plus_q_w2_w1,
      G2_DN_n1_m1_k1_k1_plus_q_w1_w1, G2_UP_n1_m1_k1_k1_plus_q_w1_w1,
      G2_DN_n2_m2_k2_plus_q_k2_w2_w2, G2_UP_n2_m2_k2_plus_q_k2_w2_w2, G4_val;

  int w_nu = parameters.get_four_point_frequency_transfer();

  scalar_type sign_div_2 = scalar_type(sign) / 2.;

  for (int w2 = 0; w2 < w_VERTEX::dmn_size(); w2++) {
    int w2_ext = w_vertex_2_w_vertex_ext(w2);
    int w2_ext_plus_w_nu = w2_ext + w_nu;
    assert(std::fabs(w_VERTEX::get_elements()[w2] - w_VERTEX_EXTENDED::get_elements()[w2_ext]) <
           1.e-6);

    for (int w1 = 0; w1 < w_VERTEX::dmn_size(); w1++) {
      int w1_ext = w_vertex_2_w_vertex_ext(w1);
      int w1_ext_plus_w_nu = w1_ext + w_nu;

      for (int k2 = 0; k2 < k_dmn_t::dmn_size(); k2++) {
        // int k2_plus_q = k_cluster_type::add(k2,q_channel);
        int k2_plus_q = q_plus_(k2);

        for (int k1 = 0; k1 < k_dmn_t::dmn_size(); k1++) {
          // int k1_plus_q = k_cluster_type::add(k1,q_channel);
          int k1_plus_q = q_plus_(k1);

          for (int m1 = 0; m1 < b::dmn_size(); m1++) {
            for (int m2 = 0; m2 < b::dmn_size(); m2++) {
              for (int n1 = 0; n1 < b::dmn_size(); n1++) {
                for (int n2 = 0; n2 < b::dmn_size(); n2++) {
                  F(n1, m2, k1, k2, w1_ext, w2_ext, G2_k_k_w_w_e_DN, G2_DN_n1_m2_k1_k2_w1_w2,
                    G2_k_k_w_w_e_UP, G2_UP_n1_m2_k1_k2_w1_w2);

                  F(n2, m1, k2_plus_q, k1_plus_q, w2_ext_plus_w_nu, w1_ext_plus_w_nu,
                    G2_k_k_w_w_e_DN, G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1, G2_k_k_w_w_e_UP,
                    G2_UP_n2_m1_k2_plus_q_k1_plus_q_w2_w1);

                  F(n1, m1, k1, k1_plus_q, w1_ext, w1_ext_plus_w_nu, G2_k_k_w_w_e_DN,
                    G2_DN_n1_m1_k1_k1_plus_q_w1_w1, G2_k_k_w_w_e_UP, G2_UP_n1_m1_k1_k1_plus_q_w1_w1);

                  F(n2, m2, k2_plus_q, k2, w2_ext_plus_w_nu, w2_ext, G2_k_k_w_w_e_DN,
                    G2_DN_n2_m2_k2_plus_q_k2_w2_w2, G2_k_k_w_w_e_UP, G2_UP_n2_m2_k2_plus_q_k2_w2_w2);

                  G4_val = -(G2_DN_n1_m2_k1_k2_w1_w2 * G2_DN_n2_m1_k2_plus_q_k1_plus_q_w2_w1 +
                             G2_UP_n1_m2_k1_k2_w1_w2 * G2_UP_n2_m1_k2_plus_q_k1_plus_q_w2_w1)

                           +
                           (G2_UP_n1_m1_k1_k1_plus_q_w1_w1 + G2_DN_n1_m1_k1_k1_plus_q_w1_w1) *
                               (G2_UP_n2_m2_k2_plus_q_k2_w2_w2 + G2_DN_n2_m2_k2_plus_q_k2_w2_w2);

                  /*
                    G4 = - (G2_k_k_w_w_e_DN(n1, m2, k1, k2, w1, w2) * G2_k_k_w_w_e_DN(n2, m1,
                    k2_plus_q, k1_plus_q, w2+w_channel, w1+w_channel)
                    + G2_k_k_w_w_e_UP(n1, m2, k1, k2, w1, w2) * G2_k_k_w_w_e_UP(n2, m1, k2_plus_q,
                    k1_plus_q, w2+w_channel, w1+w_channel))

                    + (G2_k_k_w_w_e_UP(n1, m1, k1, k1_plus_q, w1, w1+w_channel) +
                    G2_k_k_w_w_e_DN(n1, m1, k1, k1_plus_q, w1, w1+w_channel))
                    * (G2_k_k_w_w_e_UP(n2, m2, k2_plus_q, k2, w2+w_channel, w2) +
                    G2_k_k_w_w_e_DN(n2, m2, k2_plus_q, k2, w2+w_channel, w2));
                    */

                  G4(n1, n2, m1, m2, k1, k2, 0, w1, w2, 0) +=
                      std::complex<double>(sign_div_2 * G4_val);
                  // MOMS.G4_k_k_w_w(n1, n2, m1, m2, k1, k2, w1, w2) +=
                  // std::complex<double>(sign_div_2 * G4);
                }
              }
            }
          }
        }
      }
    }
  }
}

template <class parameters_type, class MOMS_type>
void accumulator_nonlocal_chi<parameters_type, MOMS_type>::accumulate_particle_particle_superconducting(
    func::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& G2_k_k_w_w_e_DN,
    func::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& G2_k_k_w_w_e_UP, scalar_type sign) {
  std::complex<scalar_type> G2_UP_n1_m1_k1_k2_w1_w2, G2_DN_n1_m1_k1_k2_w1_w2,
      G2_UP_n2_m2_q_min_k1_q_min_k2_min_w1_min_w2, G2_DN_n2_m2_q_min_k1_q_min_k2_min_w1_min_w2,
      G4_val;

  int w_nu = parameters.get_four_point_frequency_transfer();

  scalar_type sign_div_2 = sign / 2.;

  for (int w2 = 0; w2 < w_VERTEX::dmn_size(); w2++) {
    int w2_ext = w_vertex_2_w_vertex_ext(w2);
    int w_nu_min_w2 = w_nu + w_vertex_2_w_vertex_ext(min_w_vertex(w2));

    for (int w1 = 0; w1 < w_VERTEX::dmn_size(); w1++) {
      int w1_ext = w_vertex_2_w_vertex_ext(w1);
      int w_nu_min_w1 = w_nu + w_vertex_2_w_vertex_ext(min_w_vertex(w1));

      for (int k1 = 0; k1 < k_dmn_t::dmn_size(); k1++) {
        int q_minus_k1 = q_min_(k1);

        for (int k2 = 0; k2 < k_dmn_t::dmn_size(); k2++) {
          int q_minus_k2 = q_min_(k2);

          for (int n1 = 0; n1 < b::dmn_size(); n1++) {
            for (int n2 = 0; n2 < b::dmn_size(); n2++) {
              for (int m1 = 0; m1 < b::dmn_size(); m1++) {
                for (int m2 = 0; m2 < b::dmn_size(); m2++) {
                  F(n1, m1, k1, k2, w1_ext, w2_ext, G2_k_k_w_w_e_DN, G2_DN_n1_m1_k1_k2_w1_w2,
                    G2_k_k_w_w_e_UP, G2_UP_n1_m1_k1_k2_w1_w2);

                  F(n2, m2, q_minus_k1, q_minus_k2, w_nu_min_w1, w_nu_min_w2, G2_k_k_w_w_e_UP,
                    G2_UP_n2_m2_q_min_k1_q_min_k2_min_w1_min_w2, G2_k_k_w_w_e_DN,
                    G2_DN_n2_m2_q_min_k1_q_min_k2_min_w1_min_w2);

                  G4_val = (G2_UP_n1_m1_k1_k2_w1_w2 * G2_DN_n2_m2_q_min_k1_q_min_k2_min_w1_min_w2 +
                            G2_DN_n1_m1_k1_k2_w1_w2 * G2_UP_n2_m2_q_min_k1_q_min_k2_min_w1_min_w2);

                  G4(n1, n2, m1, m2, k1, k2, 0, w1, w2, 0) +=
                      std::complex<double>(sign_div_2 * G4_val);
                }
              }
            }
          }
        }
      }
    }
  }
}

}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_ACCUMULATOR_TP_ACCUMULATOR_NONLOCAL_CHI_HPP
