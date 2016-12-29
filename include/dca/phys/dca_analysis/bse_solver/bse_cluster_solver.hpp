// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class computes the vertex function \Gamma on the cluster.

#ifndef DCA_PHYS_DCA_ANALYSIS_BSE_SOLVER_BSE_CLUSTER_SOLVER_HPP
#define DCA_PHYS_DCA_ANALYSIS_BSE_SOLVER_BSE_CLUSTER_SOLVER_HPP

#include <complex>
#include <iostream>
#include <stdexcept>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/dca_step/symmetrization/diagrammatic_symmetries.hpp"
#include "dca/phys/dca_step/symmetrization/symmetrize.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/e_spin_states.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"
#include "dca/phys/vertex_measurement_type.hpp"

namespace dca {
namespace phys {
namespace analysis {
// dca::phys::analysis::

template <typename ParametersType, typename DcaDataType>
class BseClusterSolver {
public:
  using scalartype = double;

  using profiler_t = typename ParametersType::profiler_type;
  using concurrency_t = typename ParametersType::concurrency_type;

  using w = func::dmn_0<domains::frequency_domain>;
  using w_VERTEX = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT>>;
  using b = func::dmn_0<domains::electron_band_domain>;
  using DCA_k_cluster_type =
      domains::cluster_domain<double, ParametersType::lattice_type::DIMENSION, domains::CLUSTER,
                              domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>;
  using k_DCA = func::dmn_0<DCA_k_cluster_type>;

  using b_b = func::dmn_variadic<b, b>;
  using cluster_eigenvector_dmn_t = func::dmn_variadic<b, b, k_DCA, w_VERTEX>;
  using DCA_matrix_dmn_t = func::dmn_variadic<cluster_eigenvector_dmn_t, cluster_eigenvector_dmn_t>;

  BseClusterSolver(ParametersType& parameters, DcaDataType& MOMS);

  template <typename Writer>
  void write(Writer& writer);

  void compute_Gamma_cluster();

  func::function<std::complex<scalartype>, DCA_matrix_dmn_t>& get_Gamma_matrix() {
    return Gamma_matrix;
  }
  func::function<std::complex<double>, func::dmn_variadic<b_b, b_b, k_DCA, w_VERTEX>>& get_G_II_0_function() {
    return G_II_0_function;
  }

private:
  void apply_symmetries_sp();
  void apply_symmetries_tp(func::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II,
                           func::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II_0);

  void load_G_II(func::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II);
  void load_G_II_0(func::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II_0);

  void load_G_II_0_function(func::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II_0);

  void solve_BSE_on_cluster(func::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II,
                            func::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II_0);

  ParametersType& parameters;
  concurrency_t& concurrency;

  DcaDataType& MOMS;

  cluster_eigenvector_dmn_t cluster_eigenvector_dmn;

  diagrammatic_symmetries<ParametersType> diagrammatic_symmetries_obj;

  func::function<std::complex<scalartype>, DCA_matrix_dmn_t> Gamma_matrix;
  func::function<std::complex<double>, func::dmn_variadic<b_b, b_b, k_DCA, w_VERTEX>> G_II_0_function;
};

template <typename ParametersType, typename DcaDataType>
BseClusterSolver<ParametersType, DcaDataType>::BseClusterSolver(ParametersType& parameters_ref,
                                                                DcaDataType& MOMS_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),
      MOMS(MOMS_ref),

      cluster_eigenvector_dmn(),

      diagrammatic_symmetries_obj(parameters),

      Gamma_matrix("Gamma_matrix"),
      G_II_0_function("G_II_0_function") {}

template <typename ParametersType, typename DcaDataType>
template <typename Writer>
void BseClusterSolver<ParametersType, DcaDataType>::write(Writer& writer) {
  writer.execute(G_II_0_function);
}

template <typename ParametersType, typename DcaDataType>
void BseClusterSolver<ParametersType, DcaDataType>::compute_Gamma_cluster() {
  func::function<std::complex<scalartype>, DCA_matrix_dmn_t> G_II("G_II");
  func::function<std::complex<scalartype>, DCA_matrix_dmn_t> G_II_0("G_II_0");

  apply_symmetries_sp();

  load_G_II(G_II);

  load_G_II_0(G_II_0);
  load_G_II_0_function(G_II_0);

  apply_symmetries_tp(G_II, G_II_0);

  solve_BSE_on_cluster(G_II, G_II_0);
}

template <typename ParametersType, typename DcaDataType>
void BseClusterSolver<ParametersType, DcaDataType>::apply_symmetries_sp() {
  if (concurrency.id() == concurrency.last())
    std::cout << "\t" << __FUNCTION__ << "\n\n";

  profiler_t prof(__FUNCTION__, __FILE__, __LINE__);

  symmetrize::execute(MOMS.Sigma, MOMS.H_symmetry);

  symmetrize::execute(MOMS.G_k_w, MOMS.H_symmetry);
}

template <typename ParametersType, typename DcaDataType>
void BseClusterSolver<ParametersType, DcaDataType>::apply_symmetries_tp(
    func::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II,
    func::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II_0) {
  if (concurrency.id() == concurrency.last())
    std::cout << "\t" << __FUNCTION__ << "\n\n";

  if (parameters.symmetrize_Gamma()) {
    if (true) {
      if (concurrency.id() == concurrency.first())
        std::cout << "symmetrize Gamma_lattice according to the symmetry-group \n" << std::endl;

      symmetrize::execute(G_II, parameters.get_four_point_momentum_transfer());
      symmetrize::execute(G_II_0, parameters.get_four_point_momentum_transfer());
    }

    if (true) {
      if (concurrency.id() == concurrency.first())
        std::cout << "symmetrize Gamma_lattice according to diagrammatic symmetries \n"
                  << std::endl;

      diagrammatic_symmetries<ParametersType> diagrammatic_symmetries_obj(parameters);

      diagrammatic_symmetries_obj.execute(G_II);
      diagrammatic_symmetries_obj.execute(G_II_0);
    }
  }
}

template <typename ParametersType, typename DcaDataType>
void BseClusterSolver<ParametersType, DcaDataType>::load_G_II(
    func::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II) {
  if (concurrency.id() == concurrency.last())
    std::cout << "\t" << __FUNCTION__ << "\n\n";

  int* coor_1 = new int[G_II.signature()];
  int* coor_2 = new int[G_II.signature()];

  for (int i = 0; i < MOMS.G4_k_k_w_w.size(); i++) {
    MOMS.G4_k_k_w_w.linind_2_subind(i, coor_2);

    // coordinate  0 1 2 3 4 5 6 7
    // G4_k_k_w_w: b b b b k k w w
    // G_II      : b b k w b b k w

    coor_1[0] = coor_2[0];
    coor_1[1] = coor_2[1];
    coor_1[2] = coor_2[4];  // k_1
    coor_1[3] = coor_2[6];  // w_1
    coor_1[4] = coor_2[2];
    coor_1[5] = coor_2[3];
    coor_1[6] = coor_2[5];  // k_2
    coor_1[7] = coor_2[7];  // w_2

    G_II(coor_1[0], coor_1[1], coor_1[2], coor_1[3], coor_1[4], coor_1[5], coor_1[6], coor_1[7]) =
        MOMS.G4_k_k_w_w(coor_2[0], coor_2[1], coor_2[2], coor_2[3], coor_2[4], coor_2[5], coor_2[6],
                        coor_2[7]);
  }

  delete[] coor_1;
  delete[] coor_2;
}

template <typename ParametersType, typename DcaDataType>
void BseClusterSolver<ParametersType, DcaDataType>::load_G_II_0(
    func::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II_0) {
  if (concurrency.id() == concurrency.last())
    std::cout << "\t" << __FUNCTION__ << "\n\n";

  G_II_0 = 0.;

  func::dmn_variadic<k_DCA, w_VERTEX> k_w_dmn;

  int W = parameters.get_sp_fermionic_frequencies();  // TODO: Replace using w::dmn_size().
  int W_vertex = w_VERTEX::dmn_size() / 2;
  int q = parameters.get_four_point_momentum_transfer_index();

  int w_nu = parameters.get_four_point_frequency_transfer();

  int coor[2];

  for (int i = 0; i < k_w_dmn.get_size(); i++) {
    k_w_dmn.linind_2_subind(i, coor);

    int k = coor[0];
    int w_vertex = coor[1];
    int w = (coor[1] - W_vertex) + W;

    int k_plus_q = k_DCA::parameter_type::add(k, q);
    int q_minus_k = k_DCA::parameter_type::subtract(k, q);

    for (int n1 = 0; n1 < b::dmn_size(); n1++) {
      for (int n2 = 0; n2 < b::dmn_size(); n2++) {
        for (int m1 = 0; m1 < b::dmn_size(); m1++) {
          for (int m2 = 0; m2 < b::dmn_size(); m2++) {
            switch (parameters.get_four_point_type()) {
              case PARTICLE_HOLE_TRANSVERSE: {
                G_II_0(n1, n2, k, w_vertex, m1, m2, k, w_vertex) =
                    -MOMS.G_k_w(n1, e_UP, m2, e_UP, k, w) *
                    MOMS.G_k_w(n2, e_UP, m1, e_UP, k_plus_q, w + w_nu);
                break;
              }

              case PARTICLE_HOLE_MAGNETIC: {
                G_II_0(n1, n2, k, w_vertex, m1, m2, k, w_vertex) =
                    -MOMS.G_k_w(n1, e_UP, m2, e_UP, k, w) *
                    MOMS.G_k_w(n2, e_UP, m1, e_UP, k_plus_q, w + w_nu);
                break;
              }
              case PARTICLE_HOLE_CHARGE: {
                G_II_0(n1, n2, k, w_vertex, m1, m2, k, w_vertex) =
                    -2. * MOMS.G_k_w(n1, e_UP, m1, e_UP, k, w) *
                    MOMS.G_k_w(n2, e_UP, m2, e_UP, k_plus_q, w + w_nu);
                break;
              }

              case PARTICLE_PARTICLE_SUPERCONDUCTING: {
                double wn = w::get_elements()[w];
                double w_nu_min_wn = w::get_elements()[w_nu + (2 * W - 1 - w)];
                double beta = parameters.get_beta();

                if (std::fabs((w_nu * M_PI / beta - wn) - w_nu_min_wn) > 1.e-6)
                  throw std::logic_error(__FUNCTION__);

                G_II_0(n1, n2, k, w_vertex, m1, m2, k, w_vertex) =
                    MOMS.G_k_w(n1, e_UP, m1, e_UP, k, w) *
                    MOMS.G_k_w(n2, e_UP, m2, e_UP, q_minus_k, w_nu + (2 * W - 1 - w));
                break;
              }

              default:
                throw std::logic_error(__FUNCTION__);
            }
          }
        }
      }
    }
  }
}

template <typename ParametersType, typename DcaDataType>
void BseClusterSolver<ParametersType, DcaDataType>::load_G_II_0_function(
    func::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II_0) {
  if (concurrency.id() == concurrency.last())
    std::cout << "\t" << __FUNCTION__ << "\n\n";

  for (int w_ind = 0; w_ind < w_VERTEX::dmn_size(); w_ind++)
    for (int K_ind = 0; K_ind < k_DCA::dmn_size(); K_ind++)

      for (int m2 = 0; m2 < b::dmn_size(); m2++)
        for (int n2 = 0; n2 < b::dmn_size(); n2++)

          for (int m1 = 0; m1 < b::dmn_size(); m1++)
            for (int n1 = 0; n1 < b::dmn_size(); n1++)
              G_II_0_function(n1, m1, n2, m2, K_ind, w_ind) =
                  G_II_0(n1, m1, K_ind, w_ind, n2, m2, K_ind, w_ind);
}

template <typename ParametersType, typename DcaDataType>
void BseClusterSolver<ParametersType, DcaDataType>::solve_BSE_on_cluster(
    func::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II,
    func::function<std::complex<scalartype>, DCA_matrix_dmn_t>& G_II_0) {
  if (concurrency.id() == concurrency.last())
    std::cout << "\t" << __FUNCTION__ << "\n\n";

  profiler_t prof(__FUNCTION__, __FILE__, __LINE__);

  if (concurrency.id() == concurrency.last())
    std::cout << "\t" << __FUNCTION__ << std::endl << std::endl;

  scalartype renorm = 1. / (parameters.get_beta() * k_DCA::dmn_size());

  int N = cluster_eigenvector_dmn.get_size();

  dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> G4_inv(N);
  dca::linalg::Matrix<std::complex<scalartype>, dca::linalg::CPU> G4_0_inv(N);

  G_II *= renorm;
  dca::linalg::matrixop::copyArrayToMatrix(N, N, &G_II(0), N, G4_inv);
  dca::linalg::matrixop::inverse(G4_inv);

  G_II_0 *= renorm;
  dca::linalg::matrixop::copyArrayToMatrix(N, N, &G_II_0(0), N, G4_0_inv);
  dca::linalg::matrixop::inverse(G4_0_inv);

  for (int j = 0; j < N; j++)
    for (int i = 0; i < N; i++)
      Gamma_matrix(i, j) = G4_0_inv(i, j) - G4_inv(i, j);
}

}  // analysis
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_ANALYSIS_BSE_SOLVER_BSE_CLUSTER_SOLVER_HPP
