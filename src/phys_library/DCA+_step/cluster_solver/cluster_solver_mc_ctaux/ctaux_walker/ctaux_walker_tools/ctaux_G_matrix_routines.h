// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class organizes the construction of the \f$ G = N\:G^{0} \f$-matrix.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_G_MATRIX_ROUTINES_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_G_MATRIX_ROUTINES_H

#include <cassert>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>

#include "comp_library/linalg/linalg.hpp"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_structs/ctaux_auxilery_field_coefficients.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_structs/ctaux_vertex_singleton.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker/ctaux_walker_tools/ctaux_G_matrix_routines/ctaux_G_matrix_routines.hpp"

namespace DCA {
namespace QMCI {
// DCA::QMCI::

template <LIN_ALG::device_type device_t, typename parameters_type>
class G_TOOLS : public G_MATRIX_TOOLS<device_t, parameters_type> {
  typedef vertex_singleton vertex_singleton_type;

  typedef typename parameters_type::concurrency_type concurrency_type;
  typedef typename parameters_type::profiler_type profiler_t;

public:
  G_TOOLS(int id, parameters_type& parameters, CV<parameters_type>& CV_obj_ref);

  double get_Gflop();

  template <class configuration_type>
  void build_G_matrix(configuration_type& full_configuration, LIN_ALG::matrix<double, device_t>& N,
                      LIN_ALG::matrix<double, device_t>& G0, LIN_ALG::matrix<double, device_t>& G,
                      e_spin_states_type e_spin);

  double compute_G_matrix_element(int configuration_e_spin_index_i, int configuration_e_spin_index_j,
                                  LIN_ALG::matrix<double, device_t>& N,
                                  LIN_ALG::matrix<double, device_t>& G_precomputed,
                                  std::vector<vertex_singleton_type>& configuration_e_spin);

  void compute_row_on_Gamma_matrix(int row_index, LIN_ALG::vector<int, device_t>& indices,
                                   LIN_ALG::vector<double, device_t>& exp_V,
                                   LIN_ALG::matrix<double, device_t>& N,
                                   LIN_ALG::matrix<double, device_t>& G_precomputed,
                                   double* result_ptr, int incr);

  void compute_col_on_Gamma_matrix(int col_index, LIN_ALG::vector<int, device_t>& indices,
                                   LIN_ALG::vector<double, device_t>& exp_V,
                                   LIN_ALG::matrix<double, device_t>& N,
                                   LIN_ALG::matrix<double, device_t>& G_precomputed,
                                   double* result_ptr, int incr);

  void compute_G_matrix_element(LIN_ALG::vector<int, LIN_ALG::CPU>& i_index,
                                LIN_ALG::vector<int, LIN_ALG::CPU>& j_index,
                                LIN_ALG::vector<bool, LIN_ALG::CPU>& is_Bennett,
                                LIN_ALG::vector<double, LIN_ALG::CPU>& exp_Vj,
                                LIN_ALG::matrix<double, device_t>& N,
                                LIN_ALG::matrix<double, device_t>& G_precomputed,
                                double* result_ptr, int incr);

private:
  double compute_G_vertex_to_old_vertex(int configuration_e_spin_index_i,
                                        int configuration_e_spin_index_j,
                                        LIN_ALG::matrix<double, device_t>& N,
                                        std::vector<vertex_singleton_type>& configuration_e_spin);

  double compute_G_vertex_to_new_vertex(int configuration_e_spin_index_i,
                                        int configuration_e_spin_index_j,
                                        LIN_ALG::matrix<double, device_t>& G);

private:
  int thread_id;
  int stream_id;

  double GFLOP;

  parameters_type& parameters;
  concurrency_type& concurrency;

  CV<parameters_type>& CV_obj;
};

template <LIN_ALG::device_type device_t, typename parameters_type>
G_TOOLS<device_t, parameters_type>::G_TOOLS(int id, parameters_type& parameters_ref,
                                            CV<parameters_type>& CV_obj_ref)
    : G_MATRIX_TOOLS<device_t, parameters_type>(id, parameters_ref),

      thread_id(id),
      stream_id(0),

      GFLOP(0.),

      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      CV_obj(CV_obj_ref) {}

template <LIN_ALG::device_type device_t, typename parameters_type>
double G_TOOLS<device_t, parameters_type>::get_Gflop() {
  double result = GFLOP;
  GFLOP = 0;

  if (result < 0)
    std::cout << __FUNCTION__ << "\t" << result << "\n";

  return result;
}

template <LIN_ALG::device_type device_t, typename parameters_type>
template <class configuration_type>
void G_TOOLS<device_t, parameters_type>::build_G_matrix(configuration_type& full_configuration,
                                                        LIN_ALG::matrix<double, device_t>& N,
                                                        LIN_ALG::matrix<double, device_t>& G0,
                                                        LIN_ALG::matrix<double, device_t>& G,
                                                        e_spin_states_type e_spin) {
  // profiler_t profiler(concurrency, "build_G_matrix", "CT-AUX", __LINE__);

  assert(full_configuration.assert_block_form(e_spin));

  std::vector<vertex_singleton_type>& configuration_e_spin = full_configuration.get(e_spin);
  int configuration_size(configuration_e_spin.size());

  // All interaction pairs are of the same spin type, which leads to a zero configuration size for
  // one of the spin types.
  if (configuration_size == 0) {
    return;
  }

  size_t vertex_index = 0;  // # of interacting-spins
  while (vertex_index < configuration_e_spin.size() &&
         configuration_e_spin[vertex_index].get_HS_spin() != HS_ZERO)
    vertex_index++;

  {
    std::pair<int, int> size = N.get_current_size();
    size.second -= vertex_index;
    G.resize_no_copy(size);
  }

  if (N.get_current_size().first > 0 && vertex_index < configuration_e_spin.size()) {
    int m = N.get_current_size().first;
    int k = G0.get_current_size().first;
    int n = G.get_current_size().second;

    int LD_N = N.get_leading_dimension();
    int LD_G0 = G0.get_leading_dimension();
    int LD_G = G.get_leading_dimension();

#ifdef AUTOTUNING_ENABLED
    std::stringstream ss;
    ss << "GEMM_1_" << int(m / 16) * 16 + 8 << "_" << int(k / 16) * 16 + 8 << "_" << n;
    profiler_t profiler_2(concurrency, ss.str().c_str(), __FILE__, __LINE__);
#endif

    LIN_ALG::GEMM<device_t>::execute('N', 'N', m, n, k, 1., N.get_ptr(0, 0), LD_N,
                                     G0.get_ptr(0, vertex_index), LD_G0, 0., G.get_ptr(0, 0), LD_G,
                                     thread_id, stream_id);

    GFLOP += 2. * double(m) * double(n) * double(k) * (1.e-9);
  }
}

/*!
 *  \f{eqnarray}{
 *    G_{i,j} &=& (N_{ij} e^{V_j} - \delta_{ij})/(e^{V_j} -1) \\
 *    G_{i,j} &=& \sum_l N(i,l) * G^{0}(l,j)
 *  \f}
 */
template <LIN_ALG::device_type device_t, typename parameters_type>
inline double G_TOOLS<device_t, parameters_type>::compute_G_matrix_element(
    int configuration_e_spin_index_i, int configuration_e_spin_index_j,
    LIN_ALG::matrix<double, device_t>& N, LIN_ALG::matrix<double, device_t>& G_precomputed,
    std::vector<vertex_singleton_type>& configuration_e_spin) {
  int vertex_index = N.get_number_of_cols() - G_precomputed.get_number_of_cols();

  double result;

  if (configuration_e_spin_index_j < vertex_index)
    result = compute_G_vertex_to_old_vertex(configuration_e_spin_index_i,
                                            configuration_e_spin_index_j, N, configuration_e_spin);
  else
    result = compute_G_vertex_to_new_vertex(
        configuration_e_spin_index_i, configuration_e_spin_index_j - vertex_index, G_precomputed);

  return result;
}

/*!
 *  \f{eqnarray}{
 *    G_{i,j} &=& (N_{ij} e^{V_j} - \delta_{ij})/(e^{V_j} -1) \\
 *    G_{i,j} &=& \sum_l N(i,l) * G^{0}(l,j)
 *  \f}
 */
template <LIN_ALG::device_type device_t, typename parameters_type>
inline void G_TOOLS<device_t, parameters_type>::compute_G_matrix_element(
    LIN_ALG::vector<int, LIN_ALG::CPU>& i_index, LIN_ALG::vector<int, LIN_ALG::CPU>& j_index,
    LIN_ALG::vector<bool, LIN_ALG::CPU>& is_Bennett, LIN_ALG::vector<double, LIN_ALG::CPU>& exp_Vj,
    LIN_ALG::matrix<double, device_t>& N, LIN_ALG::matrix<double, device_t>& G_precomputed,
    double* result_ptr, int incr) {
  G_MATRIX_TOOLS<device_t, parameters_type>::read_G_matrix_elements(
      i_index, j_index, is_Bennett, exp_Vj, N, G_precomputed, result_ptr, incr);
}

/*!
 *  \f{eqnarray}{
 *    G_{i,j} &=& (N_{ij} e^{V_j} - \delta_{ij})/(e^{V_j} -1) \\
 *    G_{i,j} &=& \sum_l N(i,l) * G^{0}(l,j)
 *  \f}
 */
template <LIN_ALG::device_type device_t, typename parameters_type>
inline void G_TOOLS<device_t, parameters_type>::compute_row_on_Gamma_matrix(
    int row_index, LIN_ALG::vector<int, device_t>& indices,
    LIN_ALG::vector<double, device_t>& exp_V, LIN_ALG::matrix<double, device_t>& N,
    LIN_ALG::matrix<double, device_t>& G_precomputed, double* result_ptr, int incr) {
  assert(row_index > -1 && row_index < indices.size());
  G_MATRIX_TOOLS<device_t, parameters_type>::compute_row_on_Gamma_matrix(
      row_index, indices, exp_V, N, G_precomputed, result_ptr, incr);
}

/*!
 *  \f{eqnarray}{
 *    G_{i,j} &=& (N_{ij} e^{V_j} - \delta_{ij})/(e^{V_j} -1) \\
 *    G_{i,j} &=& \sum_l N(i,l) * G^{0}(l,j)
 *  \f}
 */
template <LIN_ALG::device_type device_t, typename parameters_type>
inline void G_TOOLS<device_t, parameters_type>::compute_col_on_Gamma_matrix(
    int col_index, LIN_ALG::vector<int, device_t>& indices,
    LIN_ALG::vector<double, device_t>& exp_V, LIN_ALG::matrix<double, device_t>& N,
    LIN_ALG::matrix<double, device_t>& G_precomputed, double* result_ptr, int incr) {
  assert(col_index > -1 && col_index < indices.size());
  G_MATRIX_TOOLS<device_t, parameters_type>::compute_col_on_Gamma_matrix(
      col_index, indices, exp_V, N, G_precomputed, result_ptr, incr);
}

/*!
 *  \f{eqnarray}{
 *    G_{i,j} = (N_{ij} e^{V_j} - delta_{i,j})/(e^{V_j} -1) \mbox{  (eqn 33)}
 *  \f}
 */
template <LIN_ALG::device_type device_t, typename parameters_type>
inline double G_TOOLS<device_t, parameters_type>::compute_G_vertex_to_old_vertex(
    int configuration_e_spin_index_i, int configuration_e_spin_index_j,
    LIN_ALG::matrix<double, device_t>& N, std::vector<vertex_singleton_type>& configuration_e_spin) {
  double delta = (configuration_e_spin_index_i == configuration_e_spin_index_j) ? 1. : 0.;

  vertex_singleton_type& v_j = configuration_e_spin[configuration_e_spin_index_j];

  double exp_V = CV_obj.exp_V(v_j);
  double N_ij = N(configuration_e_spin_index_i, configuration_e_spin_index_j);

  return (N_ij * exp_V - delta) / (exp_V - 1.);
}

/*!
 *  \f{eqnarray}{
 *    G_{i,j} = \sum_l N(i,l)*G0(l,j)
 *  \f}
 */
template <LIN_ALG::device_type device_t, typename parameters_type>
inline double G_TOOLS<device_t, parameters_type>::compute_G_vertex_to_new_vertex(
    int configuration_e_spin_index_i, int configuration_e_spin_index_j,
    LIN_ALG::matrix<double, device_t>& G) {
  return G(configuration_e_spin_index_i, configuration_e_spin_index_j);
}

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_G_MATRIX_ROUTINES_H
