// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class organizes the construction of the N-matrix.

/*
 *  The N-matrix can be computed directly from the Hubbard-spin configuration,
 *
 *  \f{eqnarray}{
 *     N^{-1} = e^{-\gamma \: HS_{field} \: HS_{spin}} - G^{0} * [e^{- \gamma * HS_{field} *
 * HS_{spin}}-1]
 *  \f}
 *
 *  If non-intracting spins are added to the configuration, then the N-matrix needs to be updated
 * (--> see formula 46 ),
 *
 *  \f{eqnarray}{
 *    N_{i,j} &=& N_{i,j} \mbox{ if } i<n && j<n
 *    N_{i,j} &=& \sum_{k} G_0_{i,k}*(exp(-\gamma*e_spin*HS_spin(k))-1*)*N_{k,j} \mbox{ if } i \leq
 * n && j<n
 *    N_{i,j} &=& \delta_{i,j}  \mbox{ if } j \leq n
 *  \f}
 */

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_N_MATRIX_ROUTINES_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_N_MATRIX_ROUTINES_H

#include <iostream>
#include <utility>
#include <vector>

#include "comp_library/linalg/linalg.hpp"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_structs/ctaux_auxilery_field_coefficients.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_structs/ctaux_vertex_singleton.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker/ctaux_walker_tools/ctaux_N_matrix_routines/ctaux_N_matrix_routines.hpp"

namespace DCA {
namespace QMCI {
// DCA::QMCI::

template <dca::linalg::DeviceType device_t, typename parameters_type>
class N_TOOLS : public N_MATRIX_TOOLS<device_t, parameters_type> {
  const static int MAX_VERTEX_SINGLETS = 4;

  typedef vertex_singleton vertex_singleton_type;

  typedef typename parameters_type::concurrency_type concurrency_type;
  typedef typename parameters_type::profiler_type profiler_t;

public:
  N_TOOLS(int id, parameters_type& parameters, CV<parameters_type>& CV_obj_ref);

  double get_Gflop();

  template <class configuration_type>
  void build_N_matrix(configuration_type& configuration, dca::linalg::Matrix<double, device_t>& N,
                      dca::linalg::Matrix<double, device_t>& G0, e_spin_states_type e_spin);

  template <class configuration_type>
  void update_N_matrix(configuration_type& full_configuration,
                       dca::linalg::Matrix<double, device_t>& G0,
                       dca::linalg::Matrix<double, device_t>& N, e_spin_states_type e_spin);

  template <class configuration_type>
  void rebuild_N_matrix_via_Gamma_LU(configuration_type& full_configuration,
                                     dca::linalg::Matrix<double, device_t>& N,
                                     dca::linalg::Matrix<double, device_t>& Gamma_LU,
                                     dca::linalg::Matrix<double, device_t>& G,
                                     e_spin_states_type e_spin);

  template <class configuration_type>
  void check_N_matrix(configuration_type& configuration, dca::linalg::Matrix<double, device_t>& G0,
                      dca::linalg::Matrix<double, device_t>& N,
                      dca::linalg::Matrix<double, device_t>& Gamma, e_spin_states_type e_spin);

private:
  void copy_rows_N(std::vector<int>& permutation, dca::linalg::Matrix<double, device_t>& N);

  template <class configuration_type>
  void compute_G_changed_vertices_to_all_vertex(dca::linalg::Matrix<double, device_t>& N,
                                                dca::linalg::Matrix<double, device_t>& G_precomputed,
                                                configuration_type& full_configuration,
                                                e_spin_states_type e_spin);

  void compute_d_vector(std::vector<int>& permutation, dca::linalg::Matrix<double, device_t>& N,
                        std::vector<HS_spin_states_type>& spin_values,
                        std::vector<vertex_singleton_type>& configuration_e_spin,
                        dca::linalg::Vector<double, dca::linalg::CPU>& d_inv);

  void scale_rows_N(std::vector<int>& permutation,
                    dca::linalg::Vector<double, dca::linalg::CPU>& d_inv,
                    dca::linalg::Matrix<double, device_t>& N);

  template <class configuration_type>
  static bool assert_N_matrix_format(configuration_type& full_configuration);

  template <class configuration_type>
  static bool assert_that_there_are_no_Bennett_spins(configuration_type& full_configuration);

private:
  int thread_id;
  int stream_id;

  double GFLOP;

  parameters_type& parameters;
  concurrency_type& concurrency;

  CV<parameters_type>& CV_obj;

  dca::linalg::Vector<double, dca::linalg::CPU> exp_gamma_s, one_min_exp_gamma_s, d_inv,
      exp_V_minus_one_val;

  dca::linalg::Matrix<double, device_t> G;
  dca::linalg::Matrix<double, device_t> N_new_spins;
  dca::linalg::Matrix<double, device_t> G0_times_exp_V_minus_one;
};

template <dca::linalg::DeviceType device_t, typename parameters_type>
N_TOOLS<device_t, parameters_type>::N_TOOLS(int id, parameters_type& parameters_ref,
                                            CV<parameters_type>& CV_obj_ref)
    :

      N_MATRIX_TOOLS<device_t, parameters_type>(id, parameters_ref),

      thread_id(id),
      stream_id(0),

      GFLOP(0.),

      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      CV_obj(CV_obj_ref),

      exp_gamma_s("exp_gamma_s (N_TOOLS)", MAX_VERTEX_SINGLETS * parameters.get_submatrix_size()),
      one_min_exp_gamma_s("one_min_exp_gamma_s (N_TOOLS)",
                          MAX_VERTEX_SINGLETS * parameters.get_submatrix_size()),

      d_inv("d_inv (N_TOOLS)", MAX_VERTEX_SINGLETS * parameters.get_submatrix_size()),
      exp_V_minus_one_val("exp_V_minus_one_val (N_TOOLS)",
                          MAX_VERTEX_SINGLETS * parameters.get_submatrix_size()),

      //     G                       ("G (N_TOOLS)"                       ,
      //     MAX_VERTEX_SINGLETS*parameters.get_submatrix_size()),
      //     N_new_spins             ("N_new_spins (N_TOOLS)"             ,
      //     MAX_VERTEX_SINGLETS*parameters.get_submatrix_size()),
      //     G0_times_exp_V_minus_one("G0_times_exp_V_minus_one (N_TOOLS)",
      //     MAX_VERTEX_SINGLETS*parameters.get_submatrix_size())

      G("G (N_TOOLS)", std::pair<int, int>(0, 0),
        std::pair<int, int>(parameters.get_initial_matrix_size(),
                            MAX_VERTEX_SINGLETS * parameters.get_submatrix_size())),
      N_new_spins("N_new_spins (N_TOOLS)", std::pair<int, int>(0, 0),
                  std::pair<int, int>(MAX_VERTEX_SINGLETS * parameters.get_submatrix_size(),
                                      parameters.get_initial_matrix_size())),
      G0_times_exp_V_minus_one(
          "G0_times_exp_V_minus_one (N_TOOLS)", std::pair<int, int>(0, 0),
          std::pair<int, int>(MAX_VERTEX_SINGLETS * parameters.get_submatrix_size(),
                              parameters.get_initial_matrix_size())) {}

template <dca::linalg::DeviceType device_t, typename parameters_type>
double N_TOOLS<device_t, parameters_type>::get_Gflop() {
  double result = GFLOP;
  GFLOP = 0.;

  if (result < 0)
    std::cout << __FUNCTION__ << "\t" << result << "\n";

  return result;
}

/*!
 *  The N-matrix can be computed directly from the Hubbard-spin configuration,
 *
 *  \f{eqnarray}{
 *     N^{-1} = e^{-\gamma \: HS_{field} \: HS_{spin}} - G^{0} * [e^{- \gamma * HS_{field} *
 * HS_{spin}}-1]
 *  \f}
 */
template <dca::linalg::DeviceType device_t, typename parameters_type>
template <class configuration_type>
void N_TOOLS<device_t, parameters_type>::build_N_matrix(configuration_type& configuration,
                                                        dca::linalg::Matrix<double, device_t>& N,
                                                        dca::linalg::Matrix<double, device_t>& G0,
                                                        e_spin_states_type e_spin) {
  std::vector<vertex_singleton_type>& configuration_e_spin = configuration.get(e_spin);
  int configuration_size(configuration_e_spin.size());

  // All interaction pairs are of the same spin type, which leads to a zero configuration size for
  // one of the spin types.
  if (configuration_size == 0) {
    return;
  }

  exp_gamma_s.resize(configuration_size);
  one_min_exp_gamma_s.resize(configuration_size);

  N.resizeNoCopy(configuration_size);

  for (int i = 0; i < configuration_size; ++i)
    exp_gamma_s[i] = CV_obj.exp_V(configuration_e_spin[i]);

  {  // GEMD
    for (int i = 0; i < configuration_size; ++i)
      one_min_exp_gamma_s[i] = (1. - exp_gamma_s[i]);

    LIN_ALG::GEMD<device_t>::execute(
        G0, N_MATRIX_TOOLS<device_t, parameters_type>::get_device_ptr(one_min_exp_gamma_s), N,
        thread_id, stream_id);
  }

  {
    double* exp_gamma_s_ptr = N_MATRIX_TOOLS<device_t, parameters_type>::get_device_ptr(exp_gamma_s);

    dca::linalg::blas::UseDevice<device_t>::axpy(configuration_size, 1., exp_gamma_s_ptr, 1, N.ptr(),
                                                 N.leadingDimension() + 1, thread_id, stream_id);
  }

  LIN_ALG::GEINV<device_t>::execute(N);
}

/*!
 *  If non-intracting spins are added to the configuration, then the N-matrix needs to be updated
 * (--> see formula 46 ),
 *
 *  \f{eqnarray}{
 *    N_{i,j} &=& N_{i,j} \mbox{ if } i<n && j<n \\
 *    N_{i,j} &=& \sum_{k} G_0_{i,k}*(exp(-\gamma*e_spin*HS_spin(k))-1*)*N_{k,j} \mbox{ if } i \leq
 * n && j<n \\
 *    N_{i,j} &=& \delta_{i,j}  \mbox{ if } j \leq n
 *  \f}
 */
template <dca::linalg::DeviceType device_t, typename parameters_type>
template <class configuration_type>
void N_TOOLS<device_t, parameters_type>::update_N_matrix(configuration_type& configuration,
                                                         dca::linalg::Matrix<double, device_t>& G0,
                                                         dca::linalg::Matrix<double, device_t>& N,
                                                         e_spin_states_type e_spin) {
  // profiler_t profiler(concurrency, "update_N_matrix", "CT-AUX", __LINE__, true);

  std::vector<vertex_singleton_type>& configuration_e_spin = configuration.get(e_spin);
  int configuration_size = configuration_e_spin.size();

  // All interaction pairs are of the same spin type, which leads to a zero configuration size for
  // one of the spin types.
  if (configuration_size == 0) {
    return;
  }

  int first_non_interacting_vertex_index = configuration.get_first_non_interacting_spin_index(e_spin);
  int first_shuffled_vertex_index = configuration.get_first_shuffled_spin_index(e_spin);

  assert(configuration.get_changed_spin_indices().size() == 0);
  assert(configuration.get_changed_spin_indices_e_spin(e_UP).size() == 0);
  assert(configuration.get_changed_spin_indices_e_spin(e_DN).size() == 0);
  assert(configuration.assert_block_form(e_spin));
  assert(first_shuffled_vertex_index >= first_non_interacting_vertex_index);

  if (first_non_interacting_vertex_index == configuration_size) {
    assert(configuration_size == N.size().first);
    // std::cout << __FUNCTION__
    //           << "\tconfiguration_size = " << configuration_size
    //           << "\tN.size().first = " << N.size().first << std::endl;
    return;
  }

  {
    N.resize(configuration_size);  // Move after the next if block?
  }

  {  // set columns to unity ...
    // profiler_t profiler(concurrency, "(a) set columns to unity", __FUNCTION__, __LINE__, true);

    int i = first_non_interacting_vertex_index;

    int N_r = N.nrRows();
    int N_c = N.nrCols();

    int LD = N.leadingDimension();

    assert(N_r == N_c);
    LIN_ALG::LASET<device_t>::set_zero(i, N_c - i, N.ptr(0, i), LD, thread_id, stream_id);
    LIN_ALG::LASET<device_t>::set_unity(N_r - i, N_c - i, N.ptr(i, i), LD, thread_id, stream_id);
    // dca::linalg::lapack::UseDevice<device_t>::laset("A", i, N_c - i, 0., 0., N.ptr(0, i), LD,
    //                                                thread_id, stream_id);
    // dca::linalg::lapack::UseDevice<device_t>::laset("A", N_r - i, N_c - i, 1., 0., N.ptr(i, i),
    //                                                LD, thread_id, stream_id);
  }

  if (first_shuffled_vertex_index == configuration_size || first_non_interacting_vertex_index == 0)
    return;

  {  // G0_times_exp_V_minus_one <--- \sum_{l=0}^{l<vertex_index} G0_{i,l}*exp_V_minus_one[l]
    // profiler_t profiler(concurrency, "(b) GEMD", __FUNCTION__, __LINE__, true);

    std::pair<int, int> size(configuration_size - first_shuffled_vertex_index,
                             first_non_interacting_vertex_index);

    G0_times_exp_V_minus_one.resizeNoCopy(size);

    exp_V_minus_one_val.resize(first_non_interacting_vertex_index);
    for (int j = 0; j < first_non_interacting_vertex_index; ++j)
      exp_V_minus_one_val[j] = CV_obj.exp_V(configuration_e_spin[j]) - 1.;

    double* diagonal_matrix_ptr =
        N_MATRIX_TOOLS<device_t, parameters_type>::get_device_ptr(exp_V_minus_one_val);

    LIN_ALG::GEMD<device_t>::execute(
        size, &G0.ptr()[first_shuffled_vertex_index], G0.leadingDimension(), diagonal_matrix_ptr,
        // N_MATRIX_TOOLS<device_t, parameters_type>::get_device_ptr(exp_V_minus_one_val),
        G0_times_exp_V_minus_one.ptr(), G0_times_exp_V_minus_one.leadingDimension(), thread_id,
        stream_id);
  }

  {  // G0_exp_V_minus_one * N
    // profiler_t profiler(concurrency, "(c) GEMM", __FUNCTION__, __LINE__, true);

    int m = configuration_size - first_shuffled_vertex_index;
    int k = first_non_interacting_vertex_index;
    int n = first_non_interacting_vertex_index;

    int LD_G0 = G0_times_exp_V_minus_one.leadingDimension();
    int LD_N = N.leadingDimension();

    dca::linalg::blas::UseDevice<device_t>::gemm(
        "N", "N", m, n, k, 1., G0_times_exp_V_minus_one.ptr(), LD_G0, N.ptr(), LD_N, 0.,
        &N.ptr()[first_shuffled_vertex_index], LD_N, thread_id, stream_id);

    GFLOP += 2. * double(m) * double(k) * double(n) * (1.e-9);
  }
}

template <dca::linalg::DeviceType device_t, typename parameters_type>
template <class configuration_type>
void N_TOOLS<device_t, parameters_type>::rebuild_N_matrix_via_Gamma_LU(
    configuration_type& full_configuration, dca::linalg::Matrix<double, device_t>& N,
    dca::linalg::Matrix<double, device_t>& Gamma,
    dca::linalg::Matrix<double, device_t>& G_precomputed, e_spin_states_type e_spin) {
  // profiler_t profiler(concurrency, "rebuild_N_matrix_via_Gamma_LU", "CT-AUX", __LINE__, true);

  int Gamma_size = Gamma.size().first;

  if (Gamma_size == 0)
    return;

  std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_spin);
  int configuration_size = configuration_e_spin.size();  // What happens if configuration_size = 0?

  std::vector<HS_spin_states_type>& spin_values =
      full_configuration.get_changed_spin_values_e_spin(e_spin);
  std::vector<int>& permutation = full_configuration.get_changed_spin_indices_e_spin(e_spin);

  assert(spin_values.size() == permutation.size());
  assert(Gamma_size == int(permutation.size()));
  assert(N.size().first == int(configuration_size));
  assert(assert_that_there_are_no_Bennett_spins(full_configuration));

  N_MATRIX_TOOLS<device_t, parameters_type>::set_permutation(permutation);

  {  // get the rows of N corresponding to the new spins => N_new_spins
    // profiler_t profiler(concurrency, "(a) resize N && copy rows", __FUNCTION__, __LINE__, true);

    N_new_spins.resizeNoCopy(std::pair<int, int>(Gamma_size, configuration_size));

    // copy_rows_N(permutation, N);
    N_MATRIX_TOOLS<device_t, parameters_type>::copy_rows(N, N_new_spins);
  }

  {  // get the columns of G corresponding to the new spins => G_new_spins
    // profiler_t profiler(concurrency, "(b) resize G && copy cols", __FUNCTION__, __LINE__, true);

    G.resizeNoCopy(std::pair<int, int>(configuration_size, Gamma_size));

    // if(true)
    {
      std::vector<double> exp_V(permutation.size());
      for (size_t l = 0; l < permutation.size(); ++l)
        exp_V[l] = CV_obj.exp_V(configuration_e_spin[permutation[l]]);

      N_MATRIX_TOOLS<device_t, parameters_type>::compute_G_cols(exp_V, N, G_precomputed, G);
    }
    // else
    // compute_G_changed_vertices_to_all_vertex(N, G_precomputed, full_configuration, e_spin);
  }

  {  // Gamma_LU * X = N(p_k,:) --> X = Gamma_inv_times_N_new_spins ==> (stored in N_new_spins)
    // profiler_t profiler(concurrency, "(c) LU-solve", __FUNCTION__, __LINE__, true);

    dca::linalg::matrixop::trsm('L', 'U', Gamma, N_new_spins, thread_id, stream_id);
    dca::linalg::matrixop::trsm('U', 'N', Gamma, N_new_spins, thread_id, stream_id);

    GFLOP += 2. * double(Gamma_size) * double(Gamma_size) * double(configuration_size) * (1.e-9);
  }

  {  // do N - G*Gamma_inv_times_N_new_spins --> N  || DGEMM --> work-horsegg
    // profiler_t profiler(concurrency, "(d) dgemm", __FUNCTION__, __LINE__, true);

    dca::linalg::matrixop::gemm(-1., G, N_new_spins, 1., N, thread_id, stream_id);

    GFLOP +=
        2. * double(configuration_size) * double(Gamma_size) * double(configuration_size) * (1.e-9);
  }

  {  // do N*D_i --> N ( = final N !)
    // profiler_t profiler(concurrency, "(e) rescale", __FUNCTION__, __LINE__, true);

    compute_d_vector(permutation, N, spin_values, configuration_e_spin, d_inv);

    N_MATRIX_TOOLS<device_t, parameters_type>::scale_rows(N);
  }
}

/*
  template<dca::linalg::DeviceType device_t, typename parameters_type>
  inline void N_TOOLS<device_t, parameters_type>::set_data()
  {
  std::vector<double> exp_V(permutation.size());
  std::vector<double> d_vec(permutation.size());

  {
  for(size_t l=0; l<permutation.size(); ++l)
  exp_V[l] = CV_obj.exp_V(configuration_e_spin[permutation[l]]);
  }

  {
  int                 spin_orbital, spin_orbital_paired;
  double              exp_delta_V;

  HS_field_sign       HS_field_sign;
  HS_spin_states_type old_HS_spin, new_HS_spin;

  for(int i=0; i<int(permutation.size()); ++i){

  HS_spin_states_type old_HS_spin = configuration_e_spin[permutation[i]].get_HS_spin();
  HS_spin_states_type new_HS_spin = spin_values[i];

  HS_field_sign HS_field_sign = configuration_e_spin[permutation[i]].get_HS_field();

  int spin_orbital        = configuration_e_spin[permutation[i]].get_spin_orbital();
  int spin_orbital_paired = configuration_e_spin[permutation[i]].get_paired_spin_orbital();

  if(old_HS_spin == HS_ZERO)
  {
  d_vec[i] = 1./CV_obj.exp_delta_V(spin_orbital, spin_orbital_paired, new_HS_spin, old_HS_spin,
  HS_field_sign);
  }
  else
  {
  d_inv[i] = 0;
  }
  }
  }

  N_MATRIX_TOOLS<device_t, parameters_type>::set_data(permutation, exp_V, d_vec);
  }
*/

template <dca::linalg::DeviceType device_t, typename parameters_type>
inline void N_TOOLS<device_t, parameters_type>::copy_rows_N(std::vector<int>& permutation,
                                                            dca::linalg::Matrix<double, device_t>& N) {
  // profiler_t profiler(concurrency, __FUNCTION__, __FUNCTION__, __LINE__);

  if (true) {
    for (size_t i = 0; i < permutation.size(); ++i)
      LIN_ALG::COPY<device_t>::row(N, permutation[i], N_new_spins, i);
  }
}

/*!
 *  If the vertex v_j was interacting, then
 *  \f{eqnarray}{
 *    G_{i,j} = (N_{i,j} e^{V_j} - delta_{i,j})/(e^{V_j} -1)
 *  \f}
 *
 *  If the vertex v_j was non-interacting, then
 *  \f{eqnarray}{
 *    G_{i,j} = \sum_l N_{i,l}*G0_{l,j}
 *  \f}
 */
template <dca::linalg::DeviceType device_t, typename parameters_type>
template <class configuration_type>
void N_TOOLS<device_t, parameters_type>::compute_G_changed_vertices_to_all_vertex(
    dca::linalg::Matrix<double, device_t>& N, dca::linalg::Matrix<double, device_t>& G_precomputed,
    configuration_type& full_configuration, e_spin_states_type e_spin) {
  std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_spin);

  std::vector<int>& permutation = full_configuration.get_changed_spin_indices_e_spin(e_spin);

  int configuration_size = configuration_e_spin.size();
  int Gamma_size = permutation.size();

  int N_ind = N.nrCols() - G_precomputed.nrCols();

  for (int l = 0; l < Gamma_size; ++l) {
    int j_ind = permutation[l];

    if (j_ind >= N_ind)
      LIN_ALG::COPY<device_t>::execute(configuration_size, G_precomputed.ptr(0, j_ind - N_ind), 1,
                                       G.ptr(0, l), 1);
    else {
      double exp_V = CV_obj.exp_V(configuration_e_spin[j_ind]);
      double denumerator = exp_V - 1.;

      double alpha = exp_V / denumerator;

      LIN_ALG::COPY<device_t>::execute(configuration_size, N.ptr(0, j_ind), 1, G.ptr(0, l), 1);
      LIN_ALG::SCALE<device_t>::execute(configuration_size, alpha, G.ptr(0, l), 1);

      // G(j_ind,l) -= 1./denumerator;
      LIN_ALG::MEMORY_MANAGEMENT<device_t>::add(G.ptr(j_ind, l), -1. / denumerator);
    }
  }
}

template <dca::linalg::DeviceType device_t, typename parameters_type>
inline void N_TOOLS<device_t, parameters_type>::compute_d_vector(
    std::vector<int>& permutation, dca::linalg::Matrix<double, device_t>& /*N*/,
    std::vector<HS_spin_states_type>& spin_values,
    std::vector<vertex_singleton_type>& configuration_e_spin,
    dca::linalg::Vector<double, dca::linalg::CPU>& d_inv) {
  int spin_orbital, spin_orbital_paired;
  int delta_r;
  double exp_delta_V;

  HS_field_sign HS_field_sign;
  HS_spin_states_type old_HS_spin, new_HS_spin;

  d_inv.resize(permutation.size());

  std::vector<int> d_index(0);
  // std::vector<int> N_index(0);

  for (int i = 0; i < int(permutation.size()); ++i) {
    old_HS_spin = configuration_e_spin[permutation[i]].get_HS_spin();
    new_HS_spin = spin_values[i];

    HS_field_sign = configuration_e_spin[permutation[i]].get_HS_field();

    spin_orbital = configuration_e_spin[permutation[i]].get_spin_orbital();
    spin_orbital_paired = configuration_e_spin[permutation[i]].get_paired_spin_orbital();

    delta_r = configuration_e_spin[permutation[i]].get_delta_r();

    if (old_HS_spin == HS_ZERO) {
      exp_delta_V = CV_obj.exp_delta_V(spin_orbital, spin_orbital_paired, new_HS_spin, old_HS_spin,
                                       HS_field_sign, delta_r);
      d_inv[i] = 1. / exp_delta_V;
    }
    else {
      d_inv[i] = 0;

      // d_inv[i] =
      // 1./LIN_ALG::MEMORY_MANAGEMENT<device_t>::get(N.ptr(permutation[i],permutation[i]));
      // d_inv[i] = 1./N(permutation[i],permutation[i]);

      // d_index.push_back(i);
    }
  }

  // N_MATRIX_TOOLS<device_t, parameters_type>::set_d_vector(d_index, N, d_inv);

  N_MATRIX_TOOLS<device_t, parameters_type>::set_d_vector(d_inv);
}

template <dca::linalg::DeviceType device_t, typename parameters_type>
inline void N_TOOLS<device_t, parameters_type>::scale_rows_N(
    std::vector<int>& permutation, dca::linalg::Vector<double, dca::linalg::CPU>& d_inv,
    dca::linalg::Matrix<double, device_t>& N) {
  for (size_t i = 0; i < permutation.size(); ++i)
    LIN_ALG::SCALE<device_t>::row(N, d_inv[i], permutation[i]);
}

template <dca::linalg::DeviceType device_t, typename parameters_type>
template <class configuration_type>
bool N_TOOLS<device_t, parameters_type>::assert_that_there_are_no_Bennett_spins(
    configuration_type& full_configuration) {
  {
    std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_UP);
    std::vector<int>& permutation = full_configuration.get_changed_spin_indices_e_spin(e_UP);

    for (size_t i = 0; i < permutation.size(); i++) {
      int configuration_index = configuration_e_spin[permutation[i]].get_configuration_index();
      // assert(!full_configuration[configuration_index].is_Bennett());
      if (full_configuration[configuration_index].is_Bennett())
        throw std::logic_error(__FUNCTION__);
    }
  }

  {
    std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_DN);
    std::vector<int>& permutation = full_configuration.get_changed_spin_indices_e_spin(e_DN);

    for (size_t i = 0; i < permutation.size(); i++) {
      int configuration_index = configuration_e_spin[permutation[i]].get_configuration_index();
      // assert(!full_configuration[configuration_index].is_Bennett());
      if (full_configuration[configuration_index].is_Bennett())
        throw std::logic_error(__FUNCTION__);
    }
  }

  return true;
}

template <dca::linalg::DeviceType device_t, typename parameters_type>
template <class configuration_type>
void N_TOOLS<device_t, parameters_type>::check_N_matrix(
    configuration_type& configuration, dca::linalg::Matrix<double, device_t>& N,
    dca::linalg::Matrix<double, device_t>& G0, dca::linalg::Matrix<double, device_t>& /*Gamma*/,
    e_spin_states_type e_spin) {
  dca::linalg::Matrix<double, device_t> N_correct(N.size(), N.capacity());

  std::cout.precision(4);

  build_N_matrix(configuration, N_correct, G0, e_spin);

  N.difference(N_correct);

  //     if(! N.difference(N_correct))
  //       {
  //        configuration.print();
  //        configuration.print(e_spin);

  //        N_correct.print();
  //        N.print();

  //        std::cout << "\t\t e_spin : " << e_spin << std::endl;
  //        Gamma.print();

  //        for(int i=0; i<N.size(); i++)
  //          std::cout << "\t\t" << configuration[i].get_HS_spin();
  //        std::cout << "\n";

  //        //throw std::logic_error(__FUNCTION__);

  //       }
}

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_N_MATRIX_ROUTINES_H
