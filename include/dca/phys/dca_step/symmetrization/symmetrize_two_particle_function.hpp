// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class symmetrizes two-particle Greens functions according to cluster symmetries, matsubara
// frequencies and band-index symmetries.

/*
 *  \section r_and_k cluster domain
 *
 *   For each symmetry operation \f$\mathcal{S}\f$ of the cluster-domain, we have
 *
 *   \f{eqnarray*}{
 *     G(\vec{k_1}, \vec{k_2}, \vec{q}) &=& G(\mathcal{S}(\vec{k_1}), \mathcal{S}(\vec{k_2}),
 * \mathcal{S}(\vec{q})) \\
 *   \f}
 *
 *   Hence, we have
 *
 *   \f{eqnarray*}{
 *     G(\vec{k_1}, \vec{k_2}, \vec{q}) &=& \frac{1}{N_{\mathcal{S}}} \sum_{\mathcal{S}}
 * G(\mathcal{S}(\vec{k_1}), \mathcal{S}(\vec{k_2}), \mathcal{S}(\vec{q})) \\
 *   \f}
 */

#ifndef DCA_PHYS_DCA_STEP_SYMMETRIZATION_SYMMETRIZE_TWO_PARTICLE_FUNCTION_HPP
#define DCA_PHYS_DCA_STEP_SYMMETRIZATION_SYMMETRIZE_TWO_PARTICLE_FUNCTION_HPP

#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/util/vector_operations.hpp"
#include "dca/phys/domains/cluster/cluster_operations.hpp"
#include "dca/phys/domains/cluster/cluster_symmetry.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"

namespace dca {
namespace phys {
// dca::phys::

class symmetrize_two_particle_function {
public:
  using b = func::dmn_0<domains::electron_band_domain>;

protected:
  template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
  static void execute(
      func::function<scalartype, func::dmn_variadic<b, b, b, b, k_dmn_t, k_dmn_t, w_dmn_t, w_dmn_t>>& f,
      std::vector<double> Q, bool do_diff = false);

  template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
  static void execute(
      func::function<scalartype, func::dmn_variadic<func::dmn_variadic<b, b, k_dmn_t, w_dmn_t>,
                                                    func::dmn_variadic<b, b, k_dmn_t, w_dmn_t>>>& f,
      std::vector<double> Q, bool do_diff = false);

  template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
  static void execute(
      func::function<scalartype,
                     func::dmn_variadic<func::dmn_variadic<b, b, k_dmn_t, w_dmn_t>,
                                        func::dmn_variadic<b, b, k_dmn_t, w_dmn_t>, k_dmn_t>>& f,
      bool do_diff = false);

private:
  template <typename scalartype>
  static void difference(scalartype val0, scalartype val1);

  template <typename k_dmn_t>
  struct Q_vector_is_invariant {
  public:
    static bool execute(int S_ind, std::vector<double> Q);
  };

  template <typename scalartype, typename k_dmn_t>
  static void execute(
      func::function<scalartype, func::dmn_variadic<func::dmn_variadic<b, b, k_dmn_t>,
                                                    func::dmn_variadic<b, b, k_dmn_t>>>& f,
      std::vector<double> Q, bool do_diff = false);

  template <typename scalartype, typename k_dmn_t>
  static void execute(
      func::function<scalartype, func::dmn_variadic<func::dmn_variadic<b, b, k_dmn_t>,
                                                    func::dmn_variadic<b, b, k_dmn_t>, k_dmn_t>>& f);
};

template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
void symmetrize_two_particle_function::execute(
    func::function<scalartype, func::dmn_variadic<b, b, b, b, k_dmn_t, k_dmn_t, w_dmn_t, w_dmn_t>>& f,
    std::vector<double> Q, bool do_diff) {
  {
    func::function<scalartype, func::dmn_variadic<func::dmn_variadic<b, b, k_dmn_t>,
                                                  func::dmn_variadic<b, b, k_dmn_t>>>
        G2;

    for (int w1 = 0; w1 < w_dmn_t::dmn_size(); w1++) {
      for (int w2 = 0; w2 < w_dmn_t::dmn_size(); w2++) {
        for (int l1 = 0; l1 < b::dmn_size(); l1++)
          for (int l2 = 0; l2 < b::dmn_size(); l2++)

            for (int l3 = 0; l3 < b::dmn_size(); l3++)
              for (int l4 = 0; l4 < b::dmn_size(); l4++)

                for (int k1 = 0; k1 < k_dmn_t::dmn_size(); k1++)
                  for (int k2 = 0; k2 < k_dmn_t::dmn_size(); k2++)
                    G2(l1, l2, k1, l3, l4, k2) = f(l1, l2, l3, l4, k1, k2, w1, w2);

        execute(G2, Q, do_diff);

        for (int l1 = 0; l1 < b::dmn_size(); l1++)
          for (int l2 = 0; l2 < b::dmn_size(); l2++)

            for (int l3 = 0; l3 < b::dmn_size(); l3++)
              for (int l4 = 0; l4 < b::dmn_size(); l4++)

                for (int k1 = 0; k1 < k_dmn_t::dmn_size(); k1++)
                  for (int k2 = 0; k2 < k_dmn_t::dmn_size(); k2++)
                    f(l1, l2, l3, l4, k1, k2, w1, w2) = G2(l1, l2, k1, l3, l4, k2);
      }
    }
  }
}

template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
void symmetrize_two_particle_function::execute(
    func::function<scalartype, func::dmn_variadic<func::dmn_variadic<b, b, k_dmn_t, w_dmn_t>,
                                                  func::dmn_variadic<b, b, k_dmn_t, w_dmn_t>>>& f,
    std::vector<double> Q, bool do_diff) {
  {
    func::function<scalartype, func::dmn_variadic<func::dmn_variadic<b, b, k_dmn_t>,
                                                  func::dmn_variadic<b, b, k_dmn_t>>>
        G2;

    for (int w1 = 0; w1 < w_dmn_t::dmn_size(); w1++) {
      for (int w2 = 0; w2 < w_dmn_t::dmn_size(); w2++) {
        for (int l1 = 0; l1 < b::dmn_size(); l1++)
          for (int l2 = 0; l2 < b::dmn_size(); l2++)

            for (int l3 = 0; l3 < b::dmn_size(); l3++)
              for (int l4 = 0; l4 < b::dmn_size(); l4++)

                for (int k1 = 0; k1 < k_dmn_t::dmn_size(); k1++)
                  for (int k2 = 0; k2 < k_dmn_t::dmn_size(); k2++)
                    G2(l1, l2, k1, l3, l4, k2) = f(l1, l2, k1, w1, l3, l4, k2, w2);

        execute(G2, Q, do_diff);

        for (int l1 = 0; l1 < b::dmn_size(); l1++)
          for (int l2 = 0; l2 < b::dmn_size(); l2++)

            for (int l3 = 0; l3 < b::dmn_size(); l3++)
              for (int l4 = 0; l4 < b::dmn_size(); l4++)

                for (int k1 = 0; k1 < k_dmn_t::dmn_size(); k1++)
                  for (int k2 = 0; k2 < k_dmn_t::dmn_size(); k2++)
                    f(l1, l2, k1, w1, l3, l4, k2, w2) = G2(l1, l2, k1, l3, l4, k2);
      }
    }
  }
}

template <typename scalartype, typename k_dmn_t, typename w_dmn_t>
void symmetrize_two_particle_function::execute(
    func::function<scalartype, func::dmn_variadic<func::dmn_variadic<b, b, k_dmn_t, w_dmn_t>,
                                                  func::dmn_variadic<b, b, k_dmn_t, w_dmn_t>, k_dmn_t>>& f,
    bool /*do_diff*/) {
  std::cout << __FUNCTION__ << std::endl;

  func::function<scalartype, func::dmn_variadic<func::dmn_variadic<b, b, k_dmn_t>,
                                                func::dmn_variadic<b, b, k_dmn_t>, k_dmn_t>>
      G2;

  for (int w1 = 0; w1 < w_dmn_t::dmn_size(); w1++) {
    for (int w2 = 0; w2 < w_dmn_t::dmn_size(); w2++) {
      for (int l1 = 0; l1 < b::dmn_size(); l1++)
        for (int l2 = 0; l2 < b::dmn_size(); l2++)

          for (int l3 = 0; l3 < b::dmn_size(); l3++)
            for (int l4 = 0; l4 < b::dmn_size(); l4++)

              for (int k1 = 0; k1 < k_dmn_t::dmn_size(); k1++)
                for (int k2 = 0; k2 < k_dmn_t::dmn_size(); k2++)
                  for (int k3 = 0; k3 < k_dmn_t::dmn_size(); k3++)
                    G2(l1, l2, k1, l3, l4, k2, k3) = f(l1, l2, k1, w1, l3, l4, k2, w2, k3);

      execute(G2);

      for (int l1 = 0; l1 < b::dmn_size(); l1++)
        for (int l2 = 0; l2 < b::dmn_size(); l2++)

          for (int l3 = 0; l3 < b::dmn_size(); l3++)
            for (int l4 = 0; l4 < b::dmn_size(); l4++)

              for (int k1 = 0; k1 < k_dmn_t::dmn_size(); k1++)
                for (int k2 = 0; k2 < k_dmn_t::dmn_size(); k2++)
                  for (int k3 = 0; k3 < k_dmn_t::dmn_size(); k3++)
                    f(l1, l2, k1, w1, l3, l4, k2, w2, k3) = G2(l1, l2, k1, l3, l4, k2, k3);
    }
  }
}

template <typename scalartype>
void symmetrize_two_particle_function::difference(scalartype val0, scalartype val1) {
  if (std::abs(val0 - val1) > 1.e-6) {
    throw std::logic_error(__PRETTY_FUNCTION__);
  }
}

template <typename k_dmn_t>
bool symmetrize_two_particle_function::Q_vector_is_invariant<k_dmn_t>::execute(int S_ind,
                                                                               std::vector<double> Q) {
  typedef typename k_dmn_t::parameter_type k_cluster_type;
  // typedef typename k_cluster_type::sym_super_cell_dmn_t sym_super_cell_dmn_t;
  typedef
      typename domains::cluster_symmetry<k_cluster_type>::sym_super_cell_dmn_t sym_super_cell_dmn_t;

  int DIMENSION = Q.size();

  double* O_ptr = sym_super_cell_dmn_t::get_elements()[S_ind].O;

  std::vector<double> q_rec(DIMENSION, 0.);

  for (int i = 0; i < DIMENSION; ++i)
    for (int j = 0; j < DIMENSION; ++j)
      q_rec[i] += O_ptr[i + DIMENSION * j] * Q[j];

  // q_rec = k_dmn_t::parameter_type::back_inside_cluster(q_rec);
  q_rec = domains::cluster_operations::translate_inside_cluster(
      q_rec, k_dmn_t::parameter_type::get_super_basis_vectors());

  q_rec = math::util::subtract(q_rec, Q);

  bool result;

  if (math::util::l2Norm2(q_rec) < 1.e-6)
    result = true;
  else
    result = false;

  return result;
}

template <typename scalartype, typename k_dmn_t>
void symmetrize_two_particle_function::execute(
    func::function<scalartype, func::dmn_variadic<func::dmn_variadic<b, b, k_dmn_t>,
                                                  func::dmn_variadic<b, b, k_dmn_t>>>& f,
    std::vector<double> Q, bool do_diff) {
  typedef typename k_dmn_t::parameter_type k_cluster_type;

  typedef
      typename domains::cluster_symmetry<k_cluster_type>::sym_super_cell_dmn_t sym_super_cell_dmn_t;

  static func::function<std::pair<int, int>,
                        func::dmn_variadic<func::dmn_variadic<k_dmn_t, b>, sym_super_cell_dmn_t>>&
      k_symmetry_matrix = domains::cluster_symmetry<k_cluster_type>::get_symmetry_matrix();

  static func::function<scalartype, func::dmn_variadic<func::dmn_variadic<b, b, k_dmn_t>,
                                                       func::dmn_variadic<b, b, k_dmn_t>>>
      f_new;
  f_new = scalartype(0.);

  double N_symmetries = 0;

  for (int S_ind = 0; S_ind < sym_super_cell_dmn_t::dmn_size(); ++S_ind) {
    if (Q_vector_is_invariant<k_dmn_t>::execute(S_ind, Q)) {
      N_symmetries += 1.;

      for (int l1 = 0; l1 < b::dmn_size(); l1++) {
        for (int l2 = 0; l2 < b::dmn_size(); l2++) {
          for (int l3 = 0; l3 < b::dmn_size(); l3++) {
            for (int l4 = 0; l4 < b::dmn_size(); l4++) {
              for (int k0 = 0; k0 < k_dmn_t::dmn_size(); ++k0) {
                for (int k1 = 0; k1 < k_dmn_t::dmn_size(); ++k1) {
                  int K0_new = k_symmetry_matrix(k0, l2, S_ind).first;
                  int K1_new = k_symmetry_matrix(k1, l4, S_ind).first;

                  int l1_new = k_symmetry_matrix(0, l1, S_ind).second;
                  int l2_new = k_symmetry_matrix(k0, l2, S_ind).second;
                  int l3_new = k_symmetry_matrix(0, l3, S_ind).second;
                  int l4_new = k_symmetry_matrix(k1, l4, S_ind).second;

                  f_new(l1, l2, k0, l3, l4, k1) += f(l1_new, l2_new, K0_new, l3_new, l4_new, K1_new);
                }
              }
            }
          }
        }
      }
    }
  }

  assert(N_symmetries > 0.);

  f_new /= N_symmetries;

  for (int ind = 0; ind < f.size(); ++ind) {
    if (do_diff) {
      difference(f(ind), f_new(ind));
    }

    f(ind) = f_new(ind);
  }
}

template <typename scalartype, typename k_dmn_t>
void symmetrize_two_particle_function::execute(
    func::function<scalartype, func::dmn_variadic<func::dmn_variadic<b, b, k_dmn_t>,
                                                  func::dmn_variadic<b, b, k_dmn_t>, k_dmn_t>>& f) {
  typedef typename k_dmn_t::parameter_type k_cluster_type;

  typedef typename k_cluster_type::sym_super_cell_dmn_t sym_super_cell_dmn_t;

  static func::function<std::pair<int, int>,
                        func::dmn_variadic<func::dmn_variadic<k_dmn_t, b>, sym_super_cell_dmn_t>>&
      k_symmetry_matrix = domains::cluster_symmetry<k_cluster_type>::get_symmetry_matrix();

  static func::function<scalartype, func::dmn_variadic<func::dmn_variadic<b, b, k_dmn_t>,
                                                       func::dmn_variadic<b, b, k_dmn_t>, k_dmn_t>>
      f_new;
  f_new = scalartype(0.);

  for (int S_ind = 0; S_ind < sym_super_cell_dmn_t::dmn_size(); ++S_ind) {
    for (int l1 = 0; l1 < b::dmn_size(); l1++) {
      for (int l2 = 0; l2 < b::dmn_size(); l2++) {
        for (int l3 = 0; l3 < b::dmn_size(); l3++) {
          for (int l4 = 0; l4 < b::dmn_size(); l4++) {
            for (int k0 = 0; k0 < k_dmn_t::dmn_size(); ++k0) {
              for (int k1 = 0; k1 < k_dmn_t::dmn_size(); ++k1) {
                for (int k2 = 0; k2 < k_dmn_t::dmn_size(); ++k2) {
                  int K0_new = k_symmetry_matrix(k0, 0, S_ind).first;
                  int K1_new = k_symmetry_matrix(k1, 0, S_ind).first;
                  int K2_new = k_symmetry_matrix(k2, 0, S_ind).first;

                  int l1_new = k_symmetry_matrix(0, l1, S_ind).second;
                  int l2_new = k_symmetry_matrix(k0, l2, S_ind).second;
                  int l3_new = k_symmetry_matrix(0, l3, S_ind).second;
                  int l4_new = k_symmetry_matrix(k1, l4, S_ind).second;

                  f_new(l1, l2, k0, l3, l4, k1, k2) +=
                      f(l1_new, l2_new, K0_new, l3_new, l4_new, K1_new, K2_new);
                }
              }
            }
          }
        }
      }
    }
  }

  f_new /= double(sym_super_cell_dmn_t::dmn_size());

  for (int ind = 0; ind < f.size(); ++ind)
    f(ind) = f_new(ind);
}

}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_SYMMETRIZATION_SYMMETRIZE_TWO_PARTICLE_FUNCTION_HPP
