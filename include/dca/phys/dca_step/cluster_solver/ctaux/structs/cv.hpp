// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class precomputes expensive exponents.
//
// \f{eqnarray}{
//   \exp_V(\nu, \mu, \sigma_{HS}, \sigma_{HF}, \Delta r) &=& e^{-\gamma(\nu, \mu, \Delta
//   r)\:\sigma_{HS}\:\sigma_{HF}}
// \f}

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_STRUCTS_CV_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_STRUCTS_CV_HPP

#include <cmath>
#include <utility>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/domains/hs_field_sign_domain.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/domains/hs_spin_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

template <typename parameters_type>
class CV {
public:
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using CDA = ClusterDomainAliases<parameters_type::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;

  typedef RClusterDmn r_dmn_t;

  typedef HS_spin_domain HS_spin_domain_type;
  typedef HS_field_sign_domain HS_field_sign_domain_type;

  typedef func::dmn_0<HS_spin_domain_type> HS_s;
  typedef func::dmn_0<HS_field_sign_domain_type> HS_f;

  typedef func::dmn_variadic<nu, nu, r_dmn_t> nu_nu_r_dmn_t;
  typedef func::dmn_variadic<nu, nu, HS_s, HS_f, r_dmn_t> nu_nu_HS_s_HS_f_r_dmn_t;
  typedef func::dmn_variadic<nu, nu, HS_s, HS_s, HS_f, r_dmn_t> nu_nu_HS_s_HS_s_HS_f_r_dmn_t;

public:
  CV(parameters_type& parameters);

  template <class stream_type>
  void to_JSON(stream_type& ss);

  static func::function<double, nu_nu_r_dmn_t>& get_H_interaction();

  int nu_nu_HS_s_HS_f_r_DCA_dmn_index(int spin_orbital_1, int spin_orbital_2,
                                      HS_spin_states_type HS_spin, HS_field_sign_type HS_field_sign,
                                      int site);

  /*!
   *   if U_{\nu, \mu} > 0
   *   --> 1
   *   else
   *   --> e^{- \gamma \sigma_{HS}}
   */
  template <typename vertex_singleton_t>
  double get_QMC_factor(vertex_singleton_t& v, HS_spin_states_type new_HS_spin);

  double exp_V(int linind);

  /*!
   *   \f{eqnarray}{
   *   \exp_V(\nu, \mu, \sigma_{HS}, \sigma_{HF}, \Delta r) &=& e^{-\gamma(\nu, \mu, \Delta
   * r)\:\sigma_{HS}\:\sigma_{HF}}
   *   \f}
   */
  template <typename vertex_singleton_t>
  double exp_V(vertex_singleton_t& v);

  /*!
   *   \f{eqnarray}{
   *   \exp_V(\nu, \mu, \sigma_{HS}, \sigma_{HF}, \Delta r) &=& e^{-\gamma(\nu, \mu, \Delta
   * r)\:\sigma_{HS}\:\sigma_{HF}}
   *   \f}
   */
  double exp_V(int spin_orbital_1, int spin_orbital_2, HS_spin_states_type HS_spin,
               HS_field_sign_type HS_field_sign, int site);

  double exp_delta_V(int linind);

  template <typename vertex_singleton_t>
  double exp_delta_V(vertex_singleton_t& v, HS_spin_states_type new_HS_spin);

  template <typename vertex_singleton_t>
  double exp_minus_delta_V(vertex_singleton_t& v, HS_spin_states_type new_HS_spin);

  double exp_delta_V(int spin_orbital_1, int spin_orbital_2, HS_spin_states_type HS_spin_1,
                     HS_spin_states_type HS_spin_2, HS_field_sign_type HS_field_sign, int site);

  double exp_minus_delta_V(int spin_orbital_1, int spin_orbital_2, HS_spin_states_type HS_spin_1,
                           HS_spin_states_type HS_spin_2, HS_field_sign_type HS_field_sign, int site);

  template <typename MOMS_type>
  void initialize(MOMS_type& MOMS);

private:
  void initialize_gamma();
  void initialize_exp_V();
  void initialize_exp_delta_V();

private:
  parameters_type& parameters;

  double BETA;
  double K_CT_AUX;
  double BANDS;
  double FULL_CLUSTER_SIZE;
  double CORRELATED_ORBITALS;

  nu_nu_HS_s_HS_f_r_dmn_t nu_nu_HS_s_HS_f_r_dmn;

  func::function<double, nu_nu_r_dmn_t> H_interaction;

  func::function<double, nu_nu_r_dmn_t> gamma_function;
  func::function<double, nu_nu_HS_s_HS_f_r_dmn_t> exp_V_function;

  func::function<double, nu_nu_HS_s_HS_s_HS_f_r_dmn_t> exp_delta_V_function;
};

template <typename parameters_type>
CV<parameters_type>::CV(parameters_type& parameters_ref)
    : parameters(parameters_ref),

      gamma_function("gamma_function"),
      exp_V_function("exp_V_function"),
      exp_delta_V_function("exp_delta_V_function") {}

template <typename parameters_type>
template <class stream_type>
void CV<parameters_type>::to_JSON(stream_type& ss) {
  ss << ",";

  gamma_function.to_JSON(ss);
  ss << ",";

  exp_V_function.to_JSON(ss);
  ss << ",";

  //     one__div__exp_V_function_min_one_function.to_JSON(ss);
  //     ss << ",";

  exp_delta_V_function.to_JSON(ss);
}

template <typename parameters_type>
func::function<double, typename CV<parameters_type>::nu_nu_r_dmn_t>& CV<
    parameters_type>::get_H_interaction() {
  static func::function<double, nu_nu_r_dmn_t> H;
  return H;
}

template <typename parameters_type>
inline int CV<parameters_type>::nu_nu_HS_s_HS_f_r_DCA_dmn_index(int spin_orbital_1,
                                                                int spin_orbital_2,
                                                                HS_spin_states_type HS_spin,
                                                                HS_field_sign_type HS_field_sign,
                                                                int site) {
  int HS_spin_ind = HS_spin_domain::to_coordinate(HS_spin);
  int HS_field_ind = HS_field_sign_domain::to_coordinate(HS_field_sign);

  return nu_nu_HS_s_HS_f_r_dmn(spin_orbital_1, spin_orbital_2, HS_spin_ind, HS_field_ind, site);
}

template <typename parameters_type>
inline double CV<parameters_type>::exp_V(int linind) {
  return exp_V_function(linind);
}

template <typename parameters_type>
template <typename vertex_singleton_t>
inline double CV<parameters_type>::get_QMC_factor(vertex_singleton_t& v,
                                                  HS_spin_states_type new_HS_spin) {
  std::pair<int, int>& spin_orbitals = v.get_spin_orbitals();
  int& delta_r = v.get_delta_r();

  if (H_interaction(spin_orbitals.first, spin_orbitals.second, delta_r) > 1.e-3) {
    return 1.;
  }

  if (H_interaction(spin_orbitals.first, spin_orbitals.second, delta_r) < -1.e-3) {
    HS_spin_states_type old_HS_spin = v.get_HS_spin();
    HS_field_sign_type HS_field = HS_FIELD_UP;

    return std::exp(-gamma_function(spin_orbitals.first, spin_orbitals.second, delta_r) *
                    (new_HS_spin - old_HS_spin) * HS_field);
  }

  return 0.;
}

template <typename parameters_type>
template <typename vertex_singleton_t>
inline double CV<parameters_type>::exp_V(vertex_singleton_t& v) {
  HS_spin_states_type HS_spin = v.get_HS_spin();
  HS_field_sign_type HS_field = v.get_HS_field();

  int spin_orbital = v.get_spin_orbital();
  int spin_orbital_paired = v.get_paired_spin_orbital();

  int delta_r = v.get_delta_r();

  return this->exp_V(spin_orbital, spin_orbital_paired, HS_spin, HS_field, delta_r);
}

template <typename parameters_type>
inline double CV<parameters_type>::exp_V(int spin_orbital_1, int spin_orbital_2,
                                         HS_spin_states_type HS_spin,
                                         HS_field_sign_type HS_field_sign, int site) {
  int HS_spin_ind = HS_spin_domain::to_coordinate(HS_spin);
  int HS_field_ind = HS_field_sign_domain::to_coordinate(HS_field_sign);

  return exp_V_function(spin_orbital_1, spin_orbital_2, HS_spin_ind, HS_field_ind, site);
}

template <typename parameters_type>
inline double CV<parameters_type>::exp_delta_V(int linind) {
  return exp_delta_V_function(linind);
}

template <typename parameters_type>
template <typename vertex_singleton_t>
inline double CV<parameters_type>::exp_delta_V(vertex_singleton_t& v,
                                               HS_spin_states_type new_HS_spin) {
  int spin_orbital_1 = v.get_spin_orbital();
  int spin_orbital_2 = v.get_paired_spin_orbital();

  int delta_r = v.get_delta_r();

  HS_spin_states_type old_HS_spin = v.get_HS_spin();
  HS_field_sign_type HS_field = v.get_HS_field();

  return this->exp_delta_V(spin_orbital_1, spin_orbital_2, new_HS_spin, old_HS_spin, HS_field,
                           delta_r);
}

template <typename parameters_type>
template <typename vertex_singleton_t>
inline double CV<parameters_type>::exp_minus_delta_V(vertex_singleton_t& v,
                                                     HS_spin_states_type new_HS_spin) {
  int spin_orbital_1 = v.get_spin_orbital();
  int spin_orbital_2 = v.get_paired_spin_orbital();

  int delta_r = v.get_delta_r();

  HS_spin_states_type old_HS_spin = v.get_HS_spin();
  HS_field_sign_type HS_field = v.get_HS_field();

  return this->exp_minus_delta_V(spin_orbital_1, spin_orbital_2, new_HS_spin, old_HS_spin, HS_field,
                                 delta_r);
}

template <typename parameters_type>
inline double CV<parameters_type>::exp_delta_V(int spin_orbital_1, int spin_orbital_2,
                                               HS_spin_states_type HS_spin_1,
                                               HS_spin_states_type HS_spin_2,
                                               HS_field_sign_type HS_field_sign, int site) {
  int HS_spin_1_ind = HS_spin_domain::to_coordinate(HS_spin_1);
  int HS_spin_2_ind = HS_spin_domain::to_coordinate(HS_spin_2);
  int HS_field_ind = HS_field_sign_domain::to_coordinate(HS_field_sign);

  return exp_delta_V_function(spin_orbital_1, spin_orbital_2, HS_spin_1_ind, HS_spin_2_ind,
                              HS_field_ind, site);
}

template <typename parameters_type>
inline double CV<parameters_type>::exp_minus_delta_V(int spin_orbital_1, int spin_orbital_2,
                                                     HS_spin_states_type HS_spin_1,
                                                     HS_spin_states_type HS_spin_2,
                                                     HS_field_sign_type HS_field_sign, int site) {
  int HS_spin_1_ind = HS_spin_domain::to_coordinate(HS_spin_1);
  int HS_spin_2_ind = HS_spin_domain::to_coordinate(HS_spin_2);
  int HS_field_ind = HS_field_sign_domain::to_coordinate(HS_field_sign);

  return exp_delta_V_function(spin_orbital_1, spin_orbital_2, HS_spin_2_ind, HS_spin_1_ind,
                              HS_field_ind, site);
}

template <class parameters_type>
template <typename MOMS_type>
void CV<parameters_type>::initialize(MOMS_type& MOMS) {
  BETA = parameters.get_beta();
  K_CT_AUX = parameters.get_expansion_parameter_K();
  BANDS = domains::electron_band_domain::get_size();
  FULL_CLUSTER_SIZE = r_dmn_t::dmn_size();

  H_interaction = MOMS.H_interactions;
  get_H_interaction() = MOMS.H_interactions;

  CORRELATED_ORBITALS = 0;

  for (int r_j = 0; r_j < FULL_CLUSTER_SIZE; ++r_j) {
    for (int r_i = 0; r_i < FULL_CLUSTER_SIZE; ++r_i) {
      int delta_r = r_dmn_t::parameter_type::subtract(r_j, r_i);  // delta_r = r_i - r_j

      for (int nu_j = 0; nu_j < 2 * BANDS; ++nu_j) {
        for (int nu_i = 0; nu_i < 2 * BANDS; ++nu_i) {
          if (std::abs(H_interaction(nu_i, nu_j, delta_r)) > 1.e-3) {
            ++CORRELATED_ORBITALS;
          }
        }
      }
    }
  }

  CORRELATED_ORBITALS /= 2.;

  initialize_gamma();
  initialize_exp_V();
  initialize_exp_delta_V();
}

template <class parameters_type>
void CV<parameters_type>::initialize_gamma() {
  // cout << __FUNCTION__ << endl;

  for (int nu_ind_i = 0; nu_ind_i < 2 * BANDS; nu_ind_i++) {
    for (int nu_ind_j = 0; nu_ind_j < 2 * BANDS; nu_ind_j++) {
      for (int r = 0; r < FULL_CLUSTER_SIZE; r++) {
        double U_i_j_r = std::fabs(H_interaction(nu_ind_i, nu_ind_j, r));

        double coshgamma = 1. + U_i_j_r * BETA * CORRELATED_ORBITALS / (2. * K_CT_AUX);

        gamma_function(nu_ind_i, nu_ind_j, r) = acosh(coshgamma);
      }
    }
  }
}

template <class parameters_type>
void CV<parameters_type>::initialize_exp_V() {
  // cout << __FUNCTION__ << endl;

  for (int nu_ind_i = 0; nu_ind_i < 2 * BANDS; nu_ind_i++) {
    for (int nu_ind_j = 0; nu_ind_j < 2 * BANDS; nu_ind_j++) {
      for (int HS_spin_ind = 0; HS_spin_ind < 3; HS_spin_ind++) {
        for (int HS_field_ind = 0; HS_field_ind < 2; HS_field_ind++) {
          for (int r = 0; r < FULL_CLUSTER_SIZE; r++) {
            HS_spin_states_type HS_spin = HS_spin_domain_type::get_elements()[HS_spin_ind];
            HS_field_sign_type HS_field = HS_field_sign_domain_type::get_elements()[HS_field_ind];

            // if(H_interaction(nu_ind_i, nu_ind_j, r)>1.e-3)
            {
              exp_V_function(nu_ind_i, nu_ind_j, HS_spin_ind, HS_field_ind, r) =
                  std::exp(-gamma_function(nu_ind_i, nu_ind_j, r) * HS_spin * HS_field);

              //                    if(std::fabs(exp_V_function(nu_ind_i, nu_ind_j, HS_spin_ind,
              //                    HS_field_ind, r)-1.) > 1.e-16)
              //                      one__div__exp_V_function_min_one_function(nu_ind_i, nu_ind_j,
              //                      HS_spin_ind, HS_field_ind, r)
              //                        = 1./(exp_V_function(nu_ind_i, nu_ind_j, HS_spin_ind,
              //                        HS_field_ind, r)-1.);
            }

            if (H_interaction(nu_ind_i, nu_ind_j, r) < -1.e-3) {
              //                            cout << r_dmn_t::get_r_cluster()[r][0] << "\t"
              //                                 << r_dmn_t::get_r_cluster()[r][1] << "\t"
              //                                 << H_interaction(nu_ind_i, nu_ind_j, r) << endl;

              HS_field = HS_FIELD_UP;

              exp_V_function(nu_ind_i, nu_ind_j, HS_spin_ind, HS_field_ind, r) =
                  std::exp(-gamma_function(nu_ind_i, nu_ind_j, r) * HS_spin *
                           HS_field);  // gamma=1? --> exp_V=1 ?!

              //                            if(std::fabs(exp_V_function(nu_ind_i, nu_ind_j,
              //                            HS_spin_ind, HS_field_ind, r)-1.) > 1.e-16)
              //                              one__div__exp_V_function_min_one_function(nu_ind_i,
              //                              nu_ind_j, HS_spin_ind, HS_field_ind, r)
              //                                = 1./(exp_V_function(nu_ind_i, nu_ind_j,
              //                                HS_spin_ind, 1, r)-1.);
            }
          }
        }
      }
    }
  }
}

template <class parameters_type>
void CV<parameters_type>::initialize_exp_delta_V() {
  for (int nu_ind_i = 0; nu_ind_i < 2 * BANDS; nu_ind_i++) {
    for (int nu_ind_j = 0; nu_ind_j < 2 * BANDS; nu_ind_j++) {
      for (int HS_spin_1_ind = 0; HS_spin_1_ind < 3; HS_spin_1_ind++) {
        for (int HS_spin_2_ind = 0; HS_spin_2_ind < 3; HS_spin_2_ind++) {
          for (int HS_field_ind = 0; HS_field_ind < 2; HS_field_ind++) {
            for (int r = 0; r < FULL_CLUSTER_SIZE; r++) {
              HS_spin_states_type HS_spin_1 = HS_spin_domain_type::get_elements()[HS_spin_1_ind];
              HS_spin_states_type HS_spin_2 = HS_spin_domain_type::get_elements()[HS_spin_2_ind];

              HS_field_sign_type HS_field = HS_field_sign_domain_type::get_elements()[HS_field_ind];

              // if(H_interaction(nu_ind_i, nu_ind_j, r)>1.e-3)
              {
                exp_delta_V_function(nu_ind_i, nu_ind_j, HS_spin_1_ind, HS_spin_2_ind, HS_field_ind,
                                     r) = std::exp(-gamma_function(nu_ind_i, nu_ind_j, r) *
                                                   (HS_spin_1 - HS_spin_2) * HS_field);
              }

              if (H_interaction(nu_ind_i, nu_ind_j, r) < -1.e-3) {
                HS_field = HS_FIELD_UP;

                exp_delta_V_function(nu_ind_i, nu_ind_j, HS_spin_1_ind, HS_spin_2_ind, HS_field_ind,
                                     r) = std::exp(-gamma_function(nu_ind_i, nu_ind_j, r) *
                                                   (HS_spin_1 - HS_spin_2) * HS_field);
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

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_STRUCTS_CV_HPP
