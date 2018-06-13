// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class computes all two-particle functions.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_GREENS_FUNCTIONS_TP_GREENS_FUNCTION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_GREENS_FUNCTIONS_TP_GREENS_FUNCTION_HPP

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/fermionic_overlap_matrices.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/fock_space.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/greens_functions/c_operator.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/greens_functions/integration_volume.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/greens_functions/tp_greens_function_data.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/hamiltonian.hpp"
#include "dca/phys/dca_step/cluster_solver/exact_diagonalization_advanced/hilbert_spaces/hilbert_space.hpp"
#include "dca/phys/domains/cluster/momentum_exchange_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_exchange_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/util/print_time.hpp"
#include "dca/util/plot.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ed {
// dca::phys::solver::ed::

template <typename parameters_type, typename ed_options>
class TpGreensFunction {
public:
  // typedef ED_type_definitions<parameters_type, b_dmn, s_dmn, RClusterDmn> ED_type_def;

  typedef typename ed_options::b_dmn b_dmn;
  typedef typename ed_options::s_dmn s_dmn;
  using CDA = ClusterDomainAliases<parameters_type::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using KClusterDmn = typename CDA::KClusterDmn;

  using WExchangeDmn = func::dmn_0<domains::FrequencyExchangeDomain>;
  using KExchangeDmn = func::dmn_0<domains::MomentumExchangeDomain>;

  typedef typename ed_options::profiler_t profiler_t;
  typedef typename ed_options::concurrency_type concurrency_type;

  typedef typename ed_options::scalar_type scalar_type;
  typedef typename ed_options::complex_type complex_type;

  typedef typename ed_options::vector_type vector_type;
  typedef typename ed_options::matrix_type matrix_type;
  typedef typename ed_options::int_matrix_type int_matrix_type;

  typedef typename ed_options::nu_dmn nu_dmn;
  typedef typename ed_options::b_s_r b_s_r_dmn_type;

  typedef typename ed_options::bs_dmn_type bs_dmn_type;

  typedef typename ed_options::bsr_dmn_type bsr_dmn_type;
  typedef typename ed_options::bsk_dmn_type bsk_dmn_type;

  typedef typename ed_options::nu_nu_r_dmn_type nu_nu_r_dmn_type;

  typedef func::dmn_variadic<nu_dmn, nu_dmn, nu_dmn, nu_dmn, RClusterDmn, RClusterDmn, RClusterDmn>
      nu_nu_nu_nu_r_r_r_dmn_type;

  typedef Hamiltonian<parameters_type, ed_options> fermionic_Hamiltonian_type;
  typedef fermionic_overlap_matrices<parameters_type, ed_options> fermionic_overlap_type;

  typedef Fock_space<parameters_type, ed_options> fermionic_Fock_space_type;
  typedef Hilbert_space<parameters_type, ed_options> Hilbert_space_type;

  typedef func::dmn_0<fermionic_Fock_space_type> fermionic_Fock_dmn_type;

  typedef tp_Greens_function_data<ed_options> tp_Greens_function_data_type;

  using w_VERTEX = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT>>;
  using w_VERTEX_EXTENDED = func::dmn_0<domains::vertex_frequency_domain<domains::EXTENDED>>;

public:
  TpGreensFunction(parameters_type& parameters_ref, fermionic_Hamiltonian_type& Hamiltonian_ref,
                   fermionic_overlap_type& overlap_ref);

  template <typename Writer>
  void write(Writer& writer);

  void compute_two_particle_Greens_function(bool interacting);

  void compute_particle_particle_superconducting_A(
      func::function<complex_type,
                     func::dmn_variadic<b_dmn, b_dmn, b_dmn, b_dmn, KClusterDmn, KClusterDmn,
                                        KExchangeDmn, w_VERTEX, w_VERTEX, WExchangeDmn>>& G4);
  void compute_particle_particle_superconducting_B(
      func::function<complex_type,
                     func::dmn_variadic<b_dmn, b_dmn, b_dmn, b_dmn, KClusterDmn, KClusterDmn,
                                        KExchangeDmn, w_VERTEX, w_VERTEX, WExchangeDmn>>& G4);

  void compute_two_particle_Greens_function(
      func::function<complex_type, func::dmn_variadic<w_VERTEX_EXTENDED, w_VERTEX_EXTENDED, w_VERTEX_EXTENDED,
                                                      nu_nu_nu_nu_r_r_r_dmn_type>>& G_tp_ref);

private:
  template <typename value_type>
  inline value_type Power(value_type x, int n);

  void compute_tp_Greens_function(std::vector<tp_Greens_function_data_type>& data_vec);

  /*!
   *   Here we try to compute sp-Greens-function as defined by its definition. It might be slower,
   *   but it is a good test to see that we do the right thing for the tp-Greens function
   */
  int has_nonzero_overlap(int HS_i, int HS_j, bool is_creation, int bsr_ind);

  void get_nonzero_overlap(int HS_i, int HS_j, bool is_creation, int bsr_ind, matrix_type& matrix,
                           matrix_type& tmp);

  // <C^+ C C^+ C>
  void compute_tp_permutations_ph_channel(int bsr_0, int bsr_1, int bsr_2, int bsr_3,
                                          std::vector<std::vector<c_operator>>& tp_perms);

  // <C+ C+ C C>
  void compute_tp_permutations_pp_channel(int bsr_0, int bsr_1, int bsr_2, int bsr_3,
                                          std::vector<std::vector<c_operator>>& tp_perms);

  complex_type compute_phi_slow(scalar_type E_i, scalar_type E_j, scalar_type E_k, scalar_type E_l,
                                scalar_type w1, scalar_type w2, scalar_type w3);

  void compute_tp_Greens_function_slow(int index, scalar_type E_i, scalar_type E_j, scalar_type E_k,
                                       scalar_type E_l, complex_type factor,
                                       std::vector<c_operator>& operators,
                                       tp_Greens_function_data_type& data);

  /*!
   *  new functions ...
   */
  //       void compute_real_space_Greens_functions(bool interacting);

  //       void renormalize_real_space_Greens_functions(bool interacting);

  //       void compute_Greens_functions_ac_slow(std::vector<tp_Greens_function_data_type>&
  //       data_vec);

  //       void compute_Greens_functions_ca_slow(std::vector<tp_Greens_function_data_type>&
  //       data_vec);

  //       void compute_tp_nonlocal_Greens_function(int bsr_i, int bsr_j,
  //                                                scalar_type  E_0,
  //                                                scalar_type  E_1,
  //                                                complex_type factor,
  //                                                tp_Greens_function_data_type& data);

  //       void compute_particle_particle_superconducting(func::function<complex_type,
  //       func::dmn_variadic<w_VERTEX_EXTENDED, w_VERTEX_EXTENDED, bsk_dmn_type, bsk_dmn_type> >&
  //       G_nonlocal,
  //                                                      func::function<complex_type,
  //                                                      func::dmn_variadic<b,b,b,b,KClusterDmn,KClusterDmn,w_VERTEX,w_VERTEX>
  //                                                      >&                            G4);

private:
  parameters_type& parameters;
  concurrency_type& concurrency;

  double CUT_OFF;

  fermionic_Hamiltonian_type& hamiltonian;
  fermionic_overlap_type& overlap;

  func::function<vector_type, fermionic_Fock_dmn_type>& eigen_energies;
  func::function<matrix_type, fermionic_Fock_dmn_type>& eigen_states;

  func::function<int, func::dmn_variadic<fermionic_Fock_dmn_type, fermionic_Fock_dmn_type,
                                         b_s_r_dmn_type>>& creation_set_all;
  func::function<int, func::dmn_variadic<fermionic_Fock_dmn_type, fermionic_Fock_dmn_type,
                                         b_s_r_dmn_type>>& annihilation_set_all;

  func::function<int, KClusterDmn> min_k_dmn_t;
  func::function<int, KClusterDmn> q_plus_;
  func::function<int, KClusterDmn> q_min_;

  func::function<int, w_VERTEX> min_w_vertex;
  func::function<int, w_VERTEX_EXTENDED> min_w_vertex_ext;

  func::function<int, w_VERTEX> w_vertex_2_w_vertex_ext;

  func::function<int, func::dmn_variadic<RClusterDmn, RClusterDmn>> rj_minus_ri;

  func::function<complex_type, func::dmn_variadic<w_VERTEX_EXTENDED, w_VERTEX_EXTENDED,
                                                  w_VERTEX_EXTENDED, nu_nu_nu_nu_r_r_r_dmn_type>>
      G_tp_non;
  func::function<complex_type, func::dmn_variadic<w_VERTEX_EXTENDED, w_VERTEX_EXTENDED,
                                                  w_VERTEX_EXTENDED, nu_nu_nu_nu_r_r_r_dmn_type>>
      G_tp_int;

  func::function<complex_type,
                 func::dmn_variadic<w_VERTEX_EXTENDED, w_VERTEX_EXTENDED, bsr_dmn_type, bsr_dmn_type>>
      G_non_w_w_r_r_nonlocal;
  func::function<complex_type,
                 func::dmn_variadic<w_VERTEX_EXTENDED, w_VERTEX_EXTENDED, bsk_dmn_type, bsk_dmn_type>>
      G_non_w_w_k_k_nonlocal;

  func::function<complex_type,
                 func::dmn_variadic<w_VERTEX_EXTENDED, w_VERTEX_EXTENDED, bsr_dmn_type, bsr_dmn_type>>
      G_int_w_w_r_r_nonlocal;
  func::function<complex_type,
                 func::dmn_variadic<w_VERTEX_EXTENDED, w_VERTEX_EXTENDED, bsk_dmn_type, bsk_dmn_type>>
      G_int_w_w_k_k_nonlocal;

  func::function<complex_type, func::dmn_variadic<b_dmn, b_dmn, b_dmn, b_dmn, KClusterDmn,
                                                  KClusterDmn, w_VERTEX, w_VERTEX>>
      G4_non;
  func::function<complex_type, func::dmn_variadic<b_dmn, b_dmn, b_dmn, b_dmn, KClusterDmn,
                                                  KClusterDmn, w_VERTEX, w_VERTEX>>
      G4_int;
};

template <typename parameters_type, typename ed_options>
TpGreensFunction<parameters_type, ed_options>::TpGreensFunction(
    parameters_type& parameters_ref, fermionic_Hamiltonian_type& Hamiltonian_ref,
    fermionic_overlap_type& overlap_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      CUT_OFF(parameters.get_eigenvalue_cut_off()),

      hamiltonian(Hamiltonian_ref),
      overlap(overlap_ref),

      eigen_energies(hamiltonian.get_eigen_energies()),
      eigen_states(hamiltonian.get_eigen_states()),

      creation_set_all(overlap.get_creation_set_all()),
      annihilation_set_all(overlap.get_annihilation_set_all()),

      min_k_dmn_t("min_k_dmn_t"),
      q_plus_("q_plus_"),
      q_min_("q_min_"),

      min_w_vertex(" min_w_vertex"),
      min_w_vertex_ext("min_w_vertex_ext"),

      w_vertex_2_w_vertex_ext("w_vertex_2_w_vertex_ext"),

      rj_minus_ri("rj_minus_ri"),

      G_tp_non("G_tp_non"),
      G_tp_int("G_tp_int"),

      G_non_w_w_r_r_nonlocal("G_non_w_w_r_r_nonlocal"),
      G_non_w_w_k_k_nonlocal("G_non_w_w_r_r_nonlocal"),

      G_int_w_w_r_r_nonlocal("G_int_w_w_r_r_nonlocal"),
      G_int_w_w_k_k_nonlocal("G_int_w_w_r_r_nonlocal") {
  {
    for (int ri = 0; ri < RClusterDmn::dmn_size(); ri++)
      for (int rj = 0; rj < RClusterDmn::dmn_size(); rj++)
        rj_minus_ri(ri, rj) = RClusterDmn::parameter_type::subtract(ri, rj);
  }

  {
    int q_channel = parameters.get_four_point_momentum_transfer_index();
    int k0_index = KClusterDmn::parameter_type::origin_index();

    for (int l = 0; l < KClusterDmn::parameter_type::get_size(); l++) {
      min_k_dmn_t(l) = KClusterDmn::parameter_type::subtract(l, k0_index);

      q_plus_(l) = KClusterDmn::parameter_type::add(l, q_channel);
      q_min_(l) = KClusterDmn::parameter_type::subtract(l, q_channel);
    }
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
}

template <typename parameters_type, typename ed_options>
template <typename Writer>
void TpGreensFunction<parameters_type, ed_options>::write(Writer& writer) {
  writer.open_group("fermionic-tp-Greens-function");

  writer.execute(G_tp_non);
  writer.execute(G_tp_int);

  writer.close_group();
}

/*!
  we have that w1+w2+w3+w4=0 and

  nu-wn = w1
  wn = w2
  -(wm) = w3
  -(nu-wm) = w4 = -(w1+w2+w3) =! -(nu-wm)

  so

  wn =
  wm =
  nu =

*/
template <typename parameters_type, typename ed_options>
void TpGreensFunction<parameters_type, ed_options>::compute_particle_particle_superconducting_A(
    func::function<complex_type,
                   func::dmn_variadic<b_dmn, b_dmn, b_dmn, b_dmn, KClusterDmn, KClusterDmn,
                                      KExchangeDmn, w_VERTEX, w_VERTEX, WExchangeDmn>>& G4) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\t" << __FUNCTION__ << std::endl;

  G4 = 0;

  for (int w_nu_idx = 0; w_nu_idx < WExchangeDmn::dmn_size(); ++w_nu_idx) {
    const int w_nu = WExchangeDmn::parameter_type::get_elements()[w_nu_idx];
    if (concurrency.id() == concurrency.first())
      std::cout << "\tw_nu : " << w_nu << std::endl;

    for (int b_0 = 0; b_0 < b_dmn::dmn_size(); b_0++) {
      for (int b_1 = 0; b_1 < b_dmn::dmn_size(); b_1++) {
        for (int b_2 = 0; b_2 < b_dmn::dmn_size(); b_2++) {
          for (int b_3 = 0; b_3 < b_dmn::dmn_size(); b_3++) {
            for (int r_0 = 0; r_0 < RClusterDmn::dmn_size(); r_0++) {
              for (int r_1 = 0; r_1 < RClusterDmn::dmn_size(); r_1++) {
                for (int r_2 = 0; r_2 < RClusterDmn::dmn_size(); r_2++) {
                  for (int wn = 0; wn < w_VERTEX::dmn_size(); wn++) {
                    int wn_ext = w_vertex_2_w_vertex_ext(wn);
                    assert(std::abs(w_VERTEX::get_elements()[wn] -
                                    w_VERTEX_EXTENDED::get_elements()[wn_ext]) < 1.e-6);

                    for (int wm = 0; wm < w_VERTEX::dmn_size(); wm++) {
                      int wm_ext = w_vertex_2_w_vertex_ext(wm);
                      assert(std::abs(w_VERTEX::get_elements()[wm] -
                                      w_VERTEX_EXTENDED::get_elements()[wm_ext]) < 1.e-6);

                      int w1 = w_nu + min_w_vertex_ext(wn_ext);
                      int w2 = wn_ext;
                      int w3 = min_w_vertex_ext(wm_ext);

                      //                         int w1 = wn_ext+w_nu;
                      //                         int w2 = min_w_vertex_ext(wn_ext);
                      //                         int w3 = wm_ext;

                      G4(b_0, b_1, b_2, b_3, 0, 0, 0, wn, wm, w_nu_idx) +=
                          G_tp_int(w1, w2, w3, b_0, 0, b_1, 1, b_2, 1, b_3, 0, 0, 0, 0);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    if (concurrency.id() == concurrency.first()) {
      std::cout << "\n";
      std::cout << 0.0 << "\t\t";
      for (int wn = 0; wn < w_VERTEX::dmn_size(); wn++)
        std::cout << w_VERTEX::get_elements()[wn] << "\t";
      std::cout << "\n\n";

      for (int wn = 0; wn < w_VERTEX::dmn_size(); wn++) {
        std::cout << w_VERTEX::get_elements()[wn] << "\t\t";
        for (int wm = 0; wm < w_VERTEX::dmn_size(); wm++)
          if (std::abs(real(G4(0, 0, 0, 0, 0, 0, 0, wn, wm, w_nu_idx))) < 1.e-10)
            std::cout << 0.0 << "\t";
          else
            std::cout << real(G4(0, 0, 0, 0, 0, 0, 0, wn, wm, w_nu_idx)) << "\t";
        std::cout << "\n";
      }
      std::cout << "\n";

      std::cout << "\n";
      std::cout << 0.0 << "\t\t";
      for (int wn = 0; wn < w_VERTEX::dmn_size(); wn++)
        std::cout << w_VERTEX::get_elements()[wn] << "\t";
      std::cout << "\n\n";

      for (int wn = 0; wn < w_VERTEX::dmn_size(); wn++) {
        std::cout << w_VERTEX::get_elements()[wn] << "\t\t";
        for (int wm = 0; wm < w_VERTEX::dmn_size(); wm++)
          if (std::abs(imag(G4(0, 0, 0, 0, 0, 0, 0, wn, wm, w_nu_idx))) < 1.e-10)
            std::cout << 0.0 << "\t";
          else
            std::cout << imag(G4(0, 0, 0, 0, 0, 0, 0, wn, wm, w_nu_idx)) << "\t";
        std::cout << "\n";
      }
      std::cout << "\n";
    }
  }
  {
    {
      std::vector<double> x, y;
      for (int wn = 0; wn < w_VERTEX::dmn_size(); wn++) {
        x.push_back(w_VERTEX::get_elements()[wn]);
        y.push_back(real(G4(0, 0, 0, 0, 0, 0, 0, wn, wn, 0)));
      }

      util::Plot::plotPoints(x, y);
    }

    {
      std::vector<double> x, y;
      for (int wn = 0; wn < w_VERTEX::dmn_size(); wn++) {
        x.push_back(w_VERTEX::get_elements()[wn]);
        y.push_back(real(G4(0, 0, 0, 0, 0, 0, 0, wn, w_VERTEX::dmn_size() - 1 - wn, 0)));
      }

      util::Plot::plotPoints(x, y);
    }
  }

  // assert(false);
}

/*!
  we have that w1+w2+w3+w4=0 and

  -(nu-wn) = w1
  -(wn   ) = w2
  (wm)    = w3
  (nu-wm) = w4 = -(w1+w2+w3) =! -(wm-nu)

  so

  wn =
  wm =
  nu =

*/
template <typename parameters_type, typename ed_options>
void TpGreensFunction<parameters_type, ed_options>::compute_particle_particle_superconducting_B(
    func::function<complex_type,
                   func::dmn_variadic<b_dmn, b_dmn, b_dmn, b_dmn, KClusterDmn, KClusterDmn,
                                      KExchangeDmn, w_VERTEX, w_VERTEX, WExchangeDmn>>& G4) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\n\n\t" << __FUNCTION__ << "\n\n";

  G4 = 0;

  for (int w_nu_idx = 0; w_nu_idx < WExchangeDmn::dmn_size(); ++w_nu_idx) {
    const int w_nu = WExchangeDmn::get_elements()[w_nu_idx];

    std::cout << "\n\n\t w_nu : " << w_nu << "\n";

    for (int b_0 = 0; b_0 < b_dmn::dmn_size(); b_0++) {
      for (int b_1 = 0; b_1 < b_dmn::dmn_size(); b_1++) {
        for (int b_2 = 0; b_2 < b_dmn::dmn_size(); b_2++) {
          for (int b_3 = 0; b_3 < b_dmn::dmn_size(); b_3++) {
            for (int r_0 = 0; r_0 < RClusterDmn::dmn_size(); r_0++) {
              for (int r_1 = 0; r_1 < RClusterDmn::dmn_size(); r_1++) {
                for (int r_2 = 0; r_2 < RClusterDmn::dmn_size(); r_2++) {
                  for (int wn = 0; wn < w_VERTEX::dmn_size(); wn++) {
                    int wn_ext = w_vertex_2_w_vertex_ext(wn);
                    assert(std::abs(w_VERTEX::get_elements()[wn] -
                                    w_VERTEX_EXTENDED::get_elements()[wn_ext]) < 1.e-6);

                    for (int wm = 0; wm < w_VERTEX::dmn_size(); wm++) {
                      int wm_ext = w_vertex_2_w_vertex_ext(wm);
                      assert(std::abs(w_VERTEX::get_elements()[wm] -
                                      w_VERTEX_EXTENDED::get_elements()[wm_ext]) < 1.e-6);

                      int w1 = wn_ext - w_nu;
                      int w2 = min_w_vertex_ext(wn_ext);
                      int w3 = wm_ext;

                      //                         int w1 = wn_ext+w_nu;
                      //                         int w2 = min_w_vertex_ext(wn_ext);
                      //                         int w3 = wm_ext;

                      // TODO check if ignoring the momentum is correct.
                      G4(b_0, b_1, b_2, b_3, 0, 0, 0, wn, wm, w_nu_idx) +=
                          G_tp_int(w1, w2, w3, b_0, 0, b_1, 1, b_2, 1, b_3, 0, 0, 0, 0);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    {
      std::cout << "\n";
      std::cout << 0.0 << "\t\t";
      for (int wn = 0; wn < w_VERTEX::dmn_size(); wn++)
        std::cout << w_VERTEX::get_elements()[wn] << "\t";
      std::cout << "\n\n";

      for (int wn = 0; wn < w_VERTEX::dmn_size(); wn++) {
        std::cout << w_VERTEX::get_elements()[wn] << "\t\t";
        for (int wm = 0; wm < w_VERTEX::dmn_size(); wm++)
          if (abs(real(G4(0, 0, 0, 0, 0, 0, 0, wn, wm, w_nu_idx))) < 1.e-10)
            std::cout << 0.0 << "\t";
          else
            std::cout << real(G4(0, 0, 0, 0, 0, 0, 0, wn, wm, w_nu_idx)) << "\t";
        std::cout << "\n";
      }
      std::cout << "\n";

      std::cout << "\n";
      std::cout << 0.0 << "\t\t";
      for (int wn = 0; wn < w_VERTEX::dmn_size(); wn++)
        std::cout << w_VERTEX::get_elements()[wn] << "\t";
      std::cout << "\n\n";

      for (int wn = 0; wn < w_VERTEX::dmn_size(); wn++) {
        std::cout << w_VERTEX::get_elements()[wn] << "\t\t";
        for (int wm = 0; wm < w_VERTEX::dmn_size(); wm++)
          if (abs(imag(G4(0, 0, 0, 0, 0, 0, 0, wn, wm, w_nu_idx))) < 1.e-10)
            std::cout << 0.0 << "\t";
          else
            std::cout << imag(G4(0, 0, 0, 0, 0, 0, 0, wn, wm, w_nu_idx)) << "\t";
        std::cout << "\n";
      }
      std::cout << "\n";
    }
  }
  {
    {
      std::vector<double> x, y;
      for (int wn = 0; wn < w_VERTEX::dmn_size(); wn++) {
        x.push_back(w_VERTEX::get_elements()[wn]);
        y.push_back(real(G4(0, 0, 0, 0, 0, 0, 0, wn, wn, 0)));
      }

      util::Plot::plotPoints(x, y);
    }

    {
      std::vector<double> x, y;
      for (int wn = 0; wn < w_VERTEX::dmn_size(); wn++) {
        x.push_back(w_VERTEX::get_elements()[wn]);
        y.push_back(real(G4(0, 0, 0, 0, 0, 0, wn, w_VERTEX::dmn_size() - 1 - wn)));
      }

      util::Plot::plotPoints(x, y);
    }
  }

  assert(false);
}

template <typename parameters_type, typename ed_options>
void TpGreensFunction<parameters_type, ed_options>::compute_two_particle_Greens_function(
    bool /*interacting*/) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\t" << __FUNCTION__ << std::endl;

  G_tp_int = 0;

  compute_two_particle_Greens_function(G_tp_int);

  //       if(interacting)
  //         compute_two_particle_Greens_function(G_tp_int);
  //       else
  //         compute_two_particle_Greens_function(G_tp_non);
}

template <typename parameters_type, typename ed_options>
void TpGreensFunction<parameters_type, ed_options>::compute_two_particle_Greens_function(
    func::function<complex_type, func::dmn_variadic<w_VERTEX_EXTENDED, w_VERTEX_EXTENDED, w_VERTEX_EXTENDED,
                                                    nu_nu_nu_nu_r_r_r_dmn_type>>& G_tp_ref) {
  if (concurrency.id() == concurrency.first())
    std::cout << "\t" << __FUNCTION__ << std::endl;

  G_tp_ref = 0;

  {
    int n_threads = 1;

    std::vector<tp_Greens_function_data_type> data_vec(n_threads);

    for (int l = 0; l < n_threads; l++)
      data_vec[l].initialize(parameters);

    compute_tp_Greens_function(data_vec);

    for (int l = 0; l < n_threads; l++)
      data_vec[l].sum_to(G_tp_ref);

    // concurrency.sum(G_tp_ref)
  }

  {
    scalar_type Z = hamiltonian.get_Z();

    G_tp_ref *= (1. / Z);

    //         scalar_type beta = parameters.get_beta();
    //         G_tp_ref *= (1./(beta*beta));
  }
}

template <typename parameters_type, typename ed_options>
template <typename value_type>
value_type TpGreensFunction<parameters_type, ed_options>::Power(value_type x, int n) {
  switch (n) {
    case 1:
      return x;
      break;

    case 2:
      return x * x;
      break;

    case 3:
      return x * x * x;
      break;

    case 4:
      return x * x * x * x;
      break;

    default:
      assert(false);
      return std::pow(x, n);
  }
}

template <typename parameters_type, typename ed_options>
int TpGreensFunction<parameters_type, ed_options>::has_nonzero_overlap(int HS_i, int HS_j,
                                                                       bool is_creation, int bsr_ind) {
  if (is_creation)
    return creation_set_all(HS_i, HS_j, bsr_ind);
  else
    return annihilation_set_all(HS_i, HS_j, bsr_ind);
}

template <typename parameters_type, typename ed_options>
void TpGreensFunction<parameters_type, ed_options>::get_nonzero_overlap(
    int HS_i, int HS_j, bool is_creation, int bsr_ind, matrix_type& matrix, matrix_type& tmp) {
  if (is_creation)
    overlap.compute_creation_matrix_fast(HS_i, HS_j, bsr_ind, matrix, tmp);
  else
    overlap.compute_annihilation_matrix_fast(HS_i, HS_j, bsr_ind, matrix, tmp);
}

template <typename parameters_type, typename ed_options>
void TpGreensFunction<parameters_type, ed_options>::compute_tp_permutations_ph_channel(
    int bsr_0, int bsr_1, int bsr_2, int bsr_3, std::vector<std::vector<c_operator>>& tp_perms) {
  tp_perms.resize(0);

  std::vector<c_operator> c_operators(4);

  {
    c_operators[0].index = 0;
    c_operators[1].index = 1;
    c_operators[2].index = 2;
    c_operators[3].index = 3;

    c_operators[0].bsr_ind = bsr_0;
    c_operators[1].bsr_ind = bsr_1;
    c_operators[2].bsr_ind = bsr_2;
    c_operators[3].bsr_ind = bsr_3;

    c_operators[0].creation = true;
    c_operators[1].creation = false;
    c_operators[2].creation = true;
    c_operators[3].creation = false;
  }

  {
    std::vector<int> indices(3);
    for (int l = 0; l < 3; l++)
      indices[l] = l;

    do {
      std::vector<c_operator> tmp;

      {
        for (int l = 0; l < 3; l++)
          tmp.push_back(c_operators[indices[l]]);

        tmp.push_back(c_operators[3]);
      }

      tp_perms.push_back(tmp);
    } while (std::next_permutation(indices.begin(), indices.end()));
  }
}

template <typename parameters_type, typename ed_options>
void TpGreensFunction<parameters_type, ed_options>::compute_tp_permutations_pp_channel(
    int bsr_0, int bsr_1, int bsr_2, int bsr_3, std::vector<std::vector<c_operator>>& tp_perms) {
  tp_perms.resize(0);

  std::vector<c_operator> c_operators(4);

  {
    c_operators[0].index = 0;
    c_operators[1].index = 1;
    c_operators[2].index = 2;
    c_operators[3].index = 3;

    c_operators[0].bsr_ind = bsr_0;
    c_operators[1].bsr_ind = bsr_1;
    c_operators[2].bsr_ind = bsr_2;
    c_operators[3].bsr_ind = bsr_3;

    c_operators[0].creation = true;
    c_operators[1].creation = true;
    c_operators[2].creation = false;
    c_operators[3].creation = false;
  }

  {
    std::vector<int> indices(3);
    for (int l = 0; l < 3; l++)
      indices[l] = l;

    do {
      std::vector<c_operator> tmp;

      {
        for (int l = 0; l < 3; l++)
          tmp.push_back(c_operators[indices[l]]);

        tmp.push_back(c_operators[3]);
      }

      tp_perms.push_back(tmp);
    } while (std::next_permutation(indices.begin(), indices.end()));
  }
}

/*!
 *    \int_0^{\beta} dt_1 \int_0^{t_1} dt_2 \int_0^{t_2} dt_3 e^{i*(a_1*t1+a_2*t2+a_3*t3)}
 */
template <typename parameters_type, typename ed_options>
typename TpGreensFunction<parameters_type, ed_options>::complex_type TpGreensFunction<
    parameters_type, ed_options>::compute_phi_slow(scalar_type E_i, scalar_type E_j,
                                                   scalar_type E_k, scalar_type E_l, scalar_type w1,
                                                   scalar_type w2, scalar_type w3) {
  bool do_regular = true;
  bool do_special = true;

  complex_type ONE(1, 0);
  complex_type I(0, 1);

  scalar_type beta = parameters.get_beta();

  complex_type a1 = (E_i - E_j) + w1 * I;
  complex_type a2 = (E_j - E_k) + w2 * I;
  complex_type a3 = (E_k - E_l) + w3 * I;

  complex_type result = 0;

  int index = 0;

  if (do_regular and abs(a1 + a2) > ed_options::get_epsilon() and
      abs(a2 + a3) > ed_options::get_epsilon() and abs(a1 + a2 + a3) > ed_options::get_epsilon()) {
    index++;

    result = (-(a1 * a3 * (a1 + a3) * std::exp(a1 * beta) * (-1. + std::exp(a2 * beta))) +
              Power(a2, 2) * (a3 * (-1. + std::exp(a1 * beta)) +
                              a1 * std::exp((a1 + a2) * beta) * (-1. + std::exp(a3 * beta))) +
              a2 * (Power(a3, 2) * (-1. + std::exp(a1 * beta)) -
                    2. * a1 * a3 * std::exp(a1 * beta) * (-1. + std::exp(a2 * beta)) +
                    Power(a1, 2) * std::exp((a1 + a2) * beta) * (-1. + std::exp(a3 * beta)))) /
             (a1 * a2 * (a1 + a2) * a3 * (a2 + a3) * (a1 + a2 + a3));

    //           complex_type result2 = (1./(I*w3+E_k-E_l))
    //             *(
    //               (1./(I*(w2+w3)+E_j-E_l))*(
    //               (std::exp(-beta*E_i)+std::exp(-beta*E_j))/(I*w1+E_i-E_j)
    //                                          -(std::exp(-beta*E_i)+std::exp(-beta*E_l))/(I*(w1+w2+w3)+E_i-E_l)
    //                                          )
    //               -(1./(I*w2+E_j-E_k))*(
    //               (std::exp(-beta*E_i)+std::exp(-beta*E_j))/(I*w1+E_i-E_j)
    //                                 -(std::exp(-beta*E_i)-std::exp(-beta*E_k))/(I*(w1+w2)+E_i-E_k)
    //                                     )
    //               );

    //      if(abs(result*std::exp(-beta*E_i)-result2)>1.e-6){
    //        std::cout << result << "\t" << result2 << "\n";
    //      }
    //      else
    //        std::cout << "OK\n";
  }

  // a1 + a2 = 0
  if (do_special and abs(a1 + a2) < ed_options::get_epsilon()) {
    if (abs(a2 + a3) > ed_options::get_epsilon()) {
      index++;

      result = -((a2 * Power(a3, 2) * beta + Power(a3, 2) * (-1. + std::exp(-(a2 * beta))) +
                  Power(a2, 2) * (1. + a3 * beta - std::exp(a3 * beta))) /
                 (Power(a2, 2) * Power(a3, 2) * (a2 + a3)));
    }
    // else

    if (abs(a2 + a3) < ed_options::get_epsilon()) {
      index++;

      result =
          (2. - 2. * std::exp(a3 * beta) + a3 * beta * (1. + std::exp(a3 * beta))) / Power(a3, 3);
    }

    // return result;
  }

  // a2 + a3 = 0
  if (do_special and abs(a2 + a3) < ed_options::get_epsilon()) {
    if (abs(a1 - a3) > ed_options::get_epsilon()) {
      index++;

      result = (-(a1 * Power(a3, 2) * beta * std::exp((a1 + a3) * beta)) +
                Power(a3, 2) * std::exp(a3 * beta) * (-1. + std::exp(a1 * beta)) +
                Power(a1, 2) * std::exp(a1 * beta) * (1. + (-1. + a3 * beta) * std::exp(a3 * beta))) /
               (Power(a1, 2) * (a1 - a3) * Power(a3, 2) * std::exp(a3 * beta));
    }
    // else

    //      if(abs(a1-a3)<ed_options::get_epsilon())
    //             {
    //          index++;

    //               result = (2. - 2.*std::exp(a3*beta) + a3*beta*(1. +
    //               std::exp(a3*beta)))/Power(a3,3);
    //             }

    // return result;
  }

  // a1 + a2 + a3 = 0
  if (do_special and abs(a1 + a2 + a3) < ed_options::get_epsilon()) {
    if (abs(a2 + a3) > ed_options::get_epsilon()) {
      index++;

      result = (-(Power(a2, 2) * std::exp(a2 * beta) * (-1. + std::exp(a3 * beta))) +
                a2 * a3 * std::exp(a2 * beta) * (2. + (-2. + a2 * beta) * std::exp(a3 * beta)) +
                Power(a3, 2) * (-1. + std::exp(a2 * beta) + a2 * beta * std::exp((a2 + a3) * beta))) /
               (a2 * Power(a3, 2) * Power(a2 + a3, 2) * std::exp((a2 + a3) * beta));
    }
    // else

    if (abs(a2 + a3) < ed_options::get_epsilon()) {
      index++;

      result = (2. - 2. * a3 * beta + Power(a3, 2) * Power(beta, 2) - 2. / std::exp(a3 * beta)) /
               (2. * Power(a3, 3));
    }

    // return result;
  }

  //       if(index>1)
  //         {
  //           std::cout << "\n\t" << index << "\n";
  //           std::cout << a1 << "\n";
  //           std::cout << a2 << "\n";
  //           std::cout << a3 << "\n";
  //           std::cout << "\n";
  //         }
  //       assert(index==1);
  assert(result == result);  // check for NAN

  return result;
}

template <typename parameters_type, typename ed_options>
void TpGreensFunction<parameters_type, ed_options>::compute_tp_Greens_function_slow(
    int index, scalar_type E_i, scalar_type E_j, scalar_type E_k, scalar_type E_l,
    complex_type factor, std::vector<c_operator>& operators, tp_Greens_function_data_type& data) {
  int w[3];

  scalar_type beta = parameters.get_beta();
  complex_type w_Ei = std::exp(-beta * E_i);

  for (int w3 = 0; w3 < w_VERTEX_EXTENDED::dmn_size(); w3++) {
    for (int w2 = 0; w2 < w_VERTEX_EXTENDED::dmn_size(); w2++) {
      for (int w1 = 0; w1 < w_VERTEX_EXTENDED::dmn_size(); w1++) {
        w[operators[0].index] = w1;
        w[operators[1].index] = w2;
        w[operators[2].index] = w3;

        scalar_type w1_val = w_VERTEX_EXTENDED::get_elements()[w[0]];
        scalar_type w2_val = w_VERTEX_EXTENDED::get_elements()[w[1]];
        scalar_type w3_val = w_VERTEX_EXTENDED::get_elements()[w[2]];

        complex_type phi = compute_phi_slow(E_i, E_j, E_k, E_l, w1_val, w2_val, w3_val);

        data.G_tp(w1, w2, w3, index) += factor * w_Ei * phi;
      }
    }
  }
}

// template<typename parameters_type, typename ed_options>
// typename
// TpGreensFunction<parameters_type, ed_options>::complex_type
// TpGreensFunction<parameters_type, ed_options>::compute_phi_num(scalar_type E_i,
//                                                                           scalar_type E_j,
//                                                                           scalar_type E_k,
//                                                                           scalar_type E_l,
//                                                                           scalar_type w1,
//                                                                           scalar_type w2,
//                                                                           scalar_type w3)
// {
//   complex_type ONE(1,0);
//   complex_type I  (0,1);

//   scalar_type beta = parameters.get_beta();

//   complex_type a1 = (E_i-E_j) + w1*I;
//   complex_type a2 = (E_j-E_k) + w2*I;
//   complex_type a3 = (E_k-E_l) + w3*I;

//   complex_type result=0;

//   std::vector<integration_volume> vec;

// }

// H. Hafermann et al 2009 EPL 85 27007
template <typename parameters_type, typename ed_options>
void TpGreensFunction<parameters_type, ed_options>::compute_tp_Greens_function(
    std::vector<tp_Greens_function_data_type>& data_vec) {
  // int w_nu = parameters.get_four_point_frequency_transfer();

  int origin = KClusterDmn::parameter_type::origin_index();

  std::vector<Hilbert_space_type>& Hilbert_spaces = fermionic_Fock_dmn_type::get_elements();

  // #ifdef USE_OPENMP_DIRECTIVES
  // #pragma omp parallel for
  // #endif
  for (int HS_0 = 0; HS_0 < Hilbert_spaces.size(); ++HS_0) {
    for (int HS_1 = 0; HS_1 < Hilbert_spaces.size(); ++HS_1) {
      for (int HS_2 = 0; HS_2 < Hilbert_spaces.size(); ++HS_2) {
        for (int HS_3 = 0; HS_3 < Hilbert_spaces.size(); ++HS_3) {
          int thread_id = 0;  // omp_get_thread_num();
          //               if(thread_id==0)
          //                 std::cout << "\t" << HS_0 << "\t" << HS_1 << "\t" << HS_2 << "\t" <<
          //                 HS_3 << "\t" << dca::util::print_time() << "\n";

          for (int nu_0 = 0; nu_0 < nu_dmn::dmn_size(); nu_0++) {
            for (int nu_1 = 0; nu_1 < nu_dmn::dmn_size(); nu_1++) {
              for (int nu_2 = 0; nu_2 < nu_dmn::dmn_size(); nu_2++) {
                for (int nu_3 = 0; nu_3 < nu_dmn::dmn_size(); nu_3++) {
                  for (int r_0 = 0; r_0 < RClusterDmn::dmn_size(); r_0++) {
                    for (int r_1 = 0; r_1 < RClusterDmn::dmn_size(); r_1++) {
                      for (int r_2 = 0; r_2 < RClusterDmn::dmn_size(); r_2++) {
                        // for(int r_3=0; r_3<RClusterDmn::dmn_size(); r_3++){

                        int r_3 = origin;

                        tp_Greens_function_data_type& data = data_vec[thread_id];

                        int bsr_0 = data.nu_r_dmn(nu_0, r_0);
                        int bsr_1 = data.nu_r_dmn(nu_1, r_1);
                        int bsr_2 = data.nu_r_dmn(nu_2, r_2);
                        int bsr_3 = data.nu_r_dmn(nu_3, r_3);

                        std::vector<std::vector<c_operator>> tp_perms;

                        compute_tp_permutations_pp_channel(bsr_0, bsr_1, bsr_2, bsr_3, tp_perms);

                        //                               std::cout << "\n";
                        for (int prm_ind = 0; prm_ind < tp_perms.size(); prm_ind++) {
                          scalar_type sign = -((prm_ind % 2) - 0.5) * 2.0;

                          std::vector<c_operator>& operators = tp_perms[prm_ind];

                          if (has_nonzero_overlap(HS_0, HS_1, operators[0].creation,
                                                  operators[0].bsr_ind) != -1 &&
                              has_nonzero_overlap(HS_1, HS_2, operators[1].creation,
                                                  operators[1].bsr_ind) != -1 &&
                              has_nonzero_overlap(HS_2, HS_3, operators[2].creation,
                                                  operators[2].bsr_ind) != -1 &&
                              has_nonzero_overlap(HS_3, HS_0, operators[3].creation,
                                                  operators[3].bsr_ind) != -1) {
                            bool done = false;

                            for (int l_0 = 0; l_0 < Hilbert_spaces[HS_0].size(); ++l_0) {
                              scalar_type E_0 = eigen_energies(HS_0)[l_0];

                              for (int l_1 = 0; l_1 < Hilbert_spaces[HS_1].size(); ++l_1) {
                                scalar_type E_1 = eigen_energies(HS_1)[l_1];

                                for (int l_2 = 0; l_2 < Hilbert_spaces[HS_2].size(); ++l_2) {
                                  scalar_type E_2 = eigen_energies(HS_2)[l_2];

                                  for (int l_3 = 0; l_3 < Hilbert_spaces[HS_3].size(); ++l_3) {
                                    scalar_type E_3 = eigen_energies(HS_3)[l_3];

                                    if (not done) {
                                      get_nonzero_overlap(HS_0, HS_1, operators[0].creation,
                                                          operators[0].bsr_ind, data.overlap_0,
                                                          data.tmp);
                                      get_nonzero_overlap(HS_1, HS_2, operators[1].creation,
                                                          operators[1].bsr_ind, data.overlap_1,
                                                          data.tmp);
                                      get_nonzero_overlap(HS_2, HS_3, operators[2].creation,
                                                          operators[2].bsr_ind, data.overlap_2,
                                                          data.tmp);
                                      get_nonzero_overlap(HS_3, HS_0, operators[3].creation,
                                                          operators[3].bsr_ind, data.overlap_3,
                                                          data.tmp);

                                      done = true;
                                    }

                                    complex_type factor = sign;

                                    factor *= data.overlap_0(l_0, l_1);
                                    factor *= data.overlap_1(l_1, l_2);
                                    factor *= data.overlap_2(l_2, l_3);
                                    factor *= data.overlap_3(l_3, l_0);

                                    if (abs(factor) > CUT_OFF) {
                                      // std::cout << factor << std::endl;

                                      int index = data.nu_nu_nu_nu_r_r_r_dmn(nu_0, nu_1, nu_2, nu_3,
                                                                             r_0, r_1, r_2);

                                      compute_tp_Greens_function_slow(index, E_0, E_1, E_2, E_3,
                                                                      factor, operators, data);
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

}  // ed
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_EXACT_DIAGONALIZATION_ADVANCED_GREENS_FUNCTIONS_TP_GREENS_FUNCTION_HPP
