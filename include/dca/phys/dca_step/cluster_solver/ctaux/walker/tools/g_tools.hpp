// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class organizes the construction of the \f$ G = N\:G^{0} \f$-matrix.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G_TOOLS_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G_TOOLS_HPP

#include <cassert>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>

#include "dca/linalg/linalg.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/cv.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/vertex_singleton.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/g_matrix_tools/g_matrix_tools.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

template <dca::linalg::DeviceType device_t, typename parameters_type>
class G_TOOLS : public G_MATRIX_TOOLS<device_t, parameters_type> {
  typedef vertex_singleton vertex_singleton_type;

  typedef typename parameters_type::concurrency_type concurrency_type;
  typedef typename parameters_type::profiler_type profiler_t;

public:
  G_TOOLS(int id, parameters_type& parameters, CV<parameters_type>& CV_obj_ref);

  double get_Gflop();

  template <class configuration_type>
  void build_G_matrix(configuration_type& full_configuration,
                      dca::linalg::Matrix<double, device_t>& N,
                      dca::linalg::Matrix<double, device_t>& G0,
                      dca::linalg::Matrix<double, device_t>& G, e_spin_states_type e_spin);

  double compute_G_matrix_element(int configuration_e_spin_index_i, int configuration_e_spin_index_j,
                                  dca::linalg::Matrix<double, device_t>& N,
                                  dca::linalg::Matrix<double, device_t>& G_precomputed,
                                  std::vector<vertex_singleton_type>& configuration_e_spin);

  void compute_row_on_Gamma_matrix(int row_index, dca::linalg::Vector<int, device_t>& indices,
                                   dca::linalg::Vector<double, device_t>& exp_V,
                                   dca::linalg::Matrix<double, device_t>& N,
                                   dca::linalg::Matrix<double, device_t>& G_precomputed,
                                   double* result_ptr, int incr);

  void compute_col_on_Gamma_matrix(int col_index, dca::linalg::Vector<int, device_t>& indices,
                                   dca::linalg::Vector<double, device_t>& exp_V,
                                   dca::linalg::Matrix<double, device_t>& N,
                                   dca::linalg::Matrix<double, device_t>& G_precomputed,
                                   double* result_ptr, int incr);

  void compute_G_matrix_element(dca::linalg::Vector<int, dca::linalg::CPU>& i_index,
                                dca::linalg::Vector<int, dca::linalg::CPU>& j_index,
                                dca::linalg::Vector<bool, dca::linalg::CPU>& is_Bennett,
                                dca::linalg::Vector<double, dca::linalg::CPU>& exp_Vj,
                                dca::linalg::Matrix<double, device_t>& N,
                                dca::linalg::Matrix<double, device_t>& G_precomputed,
                                double* result_ptr, int incr);



  private:
  double compute_G_vertex_to_old_vertex(int configuration_e_spin_index_i,
                                        int configuration_e_spin_index_j,
                                        dca::linalg::Matrix<double, device_t>& N,
                                        std::vector<vertex_singleton_type>& configuration_e_spin);

  double compute_G_vertex_to_new_vertex(int configuration_e_spin_index_i,
                                        int configuration_e_spin_index_j,
                                        dca::linalg::Matrix<double, device_t>& G);

private:
  int thread_id;
  int stream_id;

  double GFLOP;

  parameters_type& parameters;
  concurrency_type& concurrency;

  CV<parameters_type>& CV_obj;
};

template <dca::linalg::DeviceType device_t, typename parameters_type>
G_TOOLS<device_t, parameters_type>::G_TOOLS(int id, parameters_type& parameters_ref,
                                            CV<parameters_type>& CV_obj_ref)
    : G_MATRIX_TOOLS<device_t, parameters_type>(id, parameters_ref),

      thread_id(id),
      stream_id(0),

      GFLOP(0.),

      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      CV_obj(CV_obj_ref) {}

template <dca::linalg::DeviceType device_t, typename parameters_type>
double G_TOOLS<device_t, parameters_type>::get_Gflop() {
  double result = GFLOP;
  GFLOP = 0;

  if (result < 0)
    std::cout << __FUNCTION__ << "\t" << result << "\n";

  return result;
}

template <dca::linalg::DeviceType device_t, typename parameters_type>
template <class configuration_type>
void G_TOOLS<device_t, parameters_type>::build_G_matrix(configuration_type& full_configuration,
                                                        dca::linalg::Matrix<double, device_t>& N,
                                                        dca::linalg::Matrix<double, device_t>& G0,
                                                        dca::linalg::Matrix<double, device_t>& G,
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
    std::pair<int, int> size = N.size();
    size.second -= vertex_index;
    G.resizeNoCopy(size);
  }

  if (N.size().first > 0 && vertex_index < configuration_e_spin.size()) {
    int m = N.size().first;
    int k = G0.size().first;
    int n = G.size().second;

    int LD_N = N.leadingDimension();
    int LD_G0 = G0.leadingDimension();
    int LD_G = G.leadingDimension();

#ifdef DCA_WITH_AUTOTUNING
    std::stringstream ss;
    ss << "GEMM_1_" << int(m / 16) * 16 + 8 << "_" << int(k / 16) * 16 + 8 << "_" << n;
    profiler_t profiler_2(ss.str().c_str(), __FILE__, __LINE__);
#endif  // DCA_WITH_AUTOTUNING

    dca::linalg::blas::UseDevice<device_t>::gemm("N", "N", m, n, k, 1., N.ptr(0, 0), LD_N,
                                                 G0.ptr(0, vertex_index), LD_G0, 0., G.ptr(0, 0),
                                                 LD_G, thread_id, stream_id);

    GFLOP += 2. * double(m) * double(n) * double(k) * (1.e-9);
  }
}

/*!
 *  \f{eqnarray}{
 *    G_{i,j} &=& (N_{ij} e^{V_j} - \delta_{ij})/(e^{V_j} -1) \\
 *    G_{i,j} &=& \sum_l N(i,l) * G^{0}(l,j)
 *  \f}
 */
template <dca::linalg::DeviceType device_t, typename parameters_type>
inline double G_TOOLS<device_t, parameters_type>::compute_G_matrix_element(
    int configuration_e_spin_index_i, int configuration_e_spin_index_j,
    dca::linalg::Matrix<double, device_t>& N, dca::linalg::Matrix<double, device_t>& G_precomputed,
    std::vector<vertex_singleton_type>& configuration_e_spin) {
  int vertex_index = N.nrCols() - G_precomputed.nrCols();

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
template <dca::linalg::DeviceType device_t, typename parameters_type>
inline void G_TOOLS<device_t, parameters_type>::compute_G_matrix_element(
    dca::linalg::Vector<int, dca::linalg::CPU>& i_index,
    dca::linalg::Vector<int, dca::linalg::CPU>& j_index,
    dca::linalg::Vector<bool, dca::linalg::CPU>& is_Bennett,
    dca::linalg::Vector<double, dca::linalg::CPU>& exp_Vj, dca::linalg::Matrix<double, device_t>& N,
    dca::linalg::Matrix<double, device_t>& G_precomputed, double* result_ptr, int incr) {
  G_MATRIX_TOOLS<device_t, parameters_type>::read_G_matrix_elements(
      i_index, j_index, is_Bennett, exp_Vj, N, G_precomputed, result_ptr, incr);
}

/*!
 *  \f{eqnarray}{
 *    G_{i,j} &=& (N_{ij} e^{V_j} - \delta_{ij})/(e^{V_j} -1) \\
 *    G_{i,j} &=& \sum_l N(i,l) * G^{0}(l,j)
 *  \f}
 */
template <dca::linalg::DeviceType device_t, typename parameters_type>
inline void G_TOOLS<device_t, parameters_type>::compute_row_on_Gamma_matrix(
    int row_index, dca::linalg::Vector<int, device_t>& indices,
    dca::linalg::Vector<double, device_t>& exp_V, dca::linalg::Matrix<double, device_t>& N,
    dca::linalg::Matrix<double, device_t>& G_precomputed, double* result_ptr, int incr) {
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
template <dca::linalg::DeviceType device_t, typename parameters_type>
inline void G_TOOLS<device_t, parameters_type>::compute_col_on_Gamma_matrix(
    int col_index, dca::linalg::Vector<int, device_t>& indices,
    dca::linalg::Vector<double, device_t>& exp_V, dca::linalg::Matrix<double, device_t>& N,
    dca::linalg::Matrix<double, device_t>& G_precomputed, double* result_ptr, int incr) {
  assert(col_index > -1 && col_index < indices.size());
  G_MATRIX_TOOLS<device_t, parameters_type>::compute_col_on_Gamma_matrix(
      col_index, indices, exp_V, N, G_precomputed, result_ptr, incr);
}

/*!
 *  \f{eqnarray}{
 *    G_{i,j} = (N_{ij} e^{V_j} - delta_{i,j})/(e^{V_j} -1) \mbox{  (eqn 33)}
 *  \f}
 */
template <dca::linalg::DeviceType device_t, typename parameters_type>
inline double G_TOOLS<device_t, parameters_type>::compute_G_vertex_to_old_vertex(
    int configuration_e_spin_index_i, int configuration_e_spin_index_j,
    dca::linalg::Matrix<double, device_t>& N,
    std::vector<vertex_singleton_type>& configuration_e_spin) {
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
template <dca::linalg::DeviceType device_t, typename parameters_type>
inline double G_TOOLS<device_t, parameters_type>::compute_G_vertex_to_new_vertex(
    int configuration_e_spin_index_i, int configuration_e_spin_index_j,
    dca::linalg::Matrix<double, device_t>& G) {
  return G(configuration_e_spin_index_i, configuration_e_spin_index_j);
}

}  // ctaux
}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_TOOLS_G_TOOLS_HPP
