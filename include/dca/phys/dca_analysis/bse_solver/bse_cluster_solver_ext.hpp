// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Peter Doak (doakpw@ornl.gov)
//
// This class computes the vertex function \Gamma on the cluster.

#ifndef DCA_PHYS_DCA_ANALYSIS_BSE_SOLVER_BSE_CLUSTER_SOLVER_EXT_HPP
#define DCA_PHYS_DCA_ANALYSIS_BSE_SOLVER_BSE_CLUSTER_SOLVER_EXT_HPP

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
#include "dca/phys/four_point_type.hpp"

namespace dca {
namespace phys {
namespace analysis {
// dca::phys::analysis::

template <typename ParametersType, typename DcaDataType, typename ScalarType>
class BseClusterSolverExt {
public:
  using profiler_t = typename ParametersType::profiler_type;
  using concurrency_t = typename ParametersType::concurrency_type;
  using Lattice = typename ParametersType::lattice_type;
  using KExDmn = typename func::dmn_0<domains::MomentumExchangeDomain>;
  using WExDmn = typename func::dmn_0<domains::FrequencyExchangeDomain>;
  using w = func::dmn_0<domains::frequency_domain>;
  using WVertexDmn = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT>>;
  using b = func::dmn_0<domains::electron_band_domain>;
  using k_DCA =
      func::dmn_0<domains::cluster_domain<double, ParametersType::lattice_type::DIMENSION, domains::CLUSTER,
                                          domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;
  using cluster_eigenvector_dmn_t = func::dmn_variadic<b, b, k_DCA, WVertexDmn>;
  using DCA_matrix_dmn_t = func::dmn_variadic<cluster_eigenvector_dmn_t, cluster_eigenvector_dmn_t>;
  using TotalG4Dmn = func::dmn_variadic<DCA_matrix_dmn_t, KExDmn, WExDmn>;

  BseClusterSolverExt(ParametersType& parameters, DcaDataType& data);

  template <typename Writer>
  void write(Writer& writer);

  void compute_Gamma_cluster();

  auto& get_Gamma_cluster() /*const*/ {
    return Gamma_cluster;
  }

private:
  void apply_symmetries_sp();
  void apply_symmetries_tp(func::function<std::complex<ScalarType>, TotalG4Dmn>& G_II,
                           func::function<std::complex<ScalarType>, TotalG4Dmn>& G_II_0);

  void load_G_II(func::function<std::complex<ScalarType>, TotalG4Dmn>& G_II);
  void load_G_II_0(func::function<std::complex<ScalarType>, TotalG4Dmn>& G_II_0);

  void load_G_II_0_function(func::function<std::complex<ScalarType>, TotalG4Dmn>& G_II_0);

  void solve_BSE_on_cluster(func::function<std::complex<ScalarType>, TotalG4Dmn>& G_II,
                            func::function<std::complex<ScalarType>, TotalG4Dmn>& G_II_0);

  ParametersType& parameters;
  concurrency_t& concurrency;

  DcaDataType& data_;

  cluster_eigenvector_dmn_t cluster_eigenvector_dmn;

  diagrammatic_symmetries<ParametersType> diagrammatic_symmetries_obj;

  func::function<std::complex<ScalarType>,
                 func::dmn_variadic<cluster_eigenvector_dmn_t, cluster_eigenvector_dmn_t, KExDmn, WExDmn>>
      Gamma_cluster;
  func::function<std::complex<double>, func::dmn_variadic<b, b, b, b, k_DCA, WVertexDmn, KExDmn, WExDmn>>
      G_II_0_function;
};

template <typename ParametersType, typename DcaDataType, typename ScalarType>
BseClusterSolverExt<ParametersType, DcaDataType, ScalarType>::BseClusterSolverExt(
    ParametersType& parameters_ref, DcaDataType& data_ref)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),
      data_(data_ref),

      cluster_eigenvector_dmn(),

      diagrammatic_symmetries_obj(parameters),

      Gamma_cluster("Gamma_cluster"),
      G_II_0_function("G_II_0_function") {}

template <typename ParametersType, typename DcaDataType, typename ScalarType>
template <typename Writer>
void BseClusterSolverExt<ParametersType, DcaDataType, ScalarType>::write(Writer& writer) {
  writer.execute(Gamma_cluster);
  writer.execute(G_II_0_function);
}

template <typename ParametersType, typename DcaDataType, typename ScalarType>
void BseClusterSolverExt<ParametersType, DcaDataType, ScalarType>::compute_Gamma_cluster() {
  profiler_type prof(__FUNCTION__, "BseClusterSolverExt", __LINE__);

  func::function<std::complex<ScalarType>, TotalG4Dmn> G_II("G_II");
  func::function<std::complex<ScalarType>, TotalG4Dmn> G_II_0("G_II_0");

  apply_symmetries_sp();

  load_G_II(G_II);

  load_G_II_0(G_II_0);
  load_G_II_0_function(G_II_0);

  apply_symmetries_tp(G_II, G_II_0);

  solve_BSE_on_cluster(G_II, G_II_0);
}

template <typename ParametersType, typename DcaDataType, typename ScalarType>
void BseClusterSolverExt<ParametersType, DcaDataType, ScalarType>::apply_symmetries_sp() {
  profiler_type prof(__FUNCTION__, "BseClusterSolverExt", __LINE__);

  if (concurrency.id() == concurrency.first())
    std::cout << "\t" << __FUNCTION__ << "\n\n";

  profiler_t prof(__FUNCTION__, __FILE__, __LINE__);

  symmetrize::execute<Lattice>(data_.Sigma, data_.H_symmetry);
  symmetrize::execute<Lattice>(data_.G_k_w, data_.H_symmetry);
}

template <typename ParametersType, typename DcaDataType, typename ScalarType>
void BseClusterSolverExt<ParametersType, DcaDataType, ScalarType>::apply_symmetries_tp(
    func::function<std::complex<ScalarType>, TotalG4Dmn>& G_II,
    func::function<std::complex<ScalarType>, TotalG4Dmn>& G_II_0) {
  profiler_type prof(__FUNCTION__, "BseClusterSolverExt", __LINE__);

  if (concurrency.id() == concurrency.first())
    std::cout << "\t" << __FUNCTION__ << "\n\n";

  if (parameters.symmetrize_Gamma()) {
    std::vector<int> subind(3);
    if (concurrency.id() == concurrency.first())
       std::cout << "symmetrize Gamma_lattice according to the symmetry-group \n" << std::endl;
    for (int wex_ind = 0; wex_ind < WExDmn::dmn_size(); wex_ind++)
      for (int kex_ind = 0; kex_ind < KExDmn::dmn_size(); kex_ind++) {
        func::function<std::complex<ScalarType>, DCA_matrix_dmn_t> G_II_indi;
        func::function<std::complex<ScalarType>, DCA_matrix_dmn_t> G_II_0_indi;
        subind = {0, kex_ind, wex_ind};
        G_II.slice(0, subind, static_cast<std::complex<ScalarType>*>(G_II_indi.values()));
        G_II_0.slice(0, subind, static_cast<std::complex<ScalarType>*>(G_II_0_indi.values()));

        symmetrize::execute(G_II_indi, parameters.get_four_point_momentum_transfer());
        symmetrize::execute(G_II_0_indi, parameters.get_four_point_momentum_transfer());
        G_II.distribute(0, subind, static_cast<std::complex<ScalarType>*>(G_II_indi.values()));
        G_II_0.distribute(0, subind, static_cast<std::complex<ScalarType>*>(G_II_0_indi.values()));
      }
  }
  if (true) {
    if (concurrency.id() == concurrency.first())
      std::cout << "symmetrize Gamma_lattice according to diagrammatic symmetries \n" << std::endl;

    diagrammatic_symmetries<ParametersType> diagrammatic_symmetries_obj(parameters);

    diagrammatic_symmetries_obj.execute(G_II);
    diagrammatic_symmetries_obj.execute(G_II_0);
  }
}

template <typename ParametersType, typename DcaDataType, typename ScalarType>
void BseClusterSolverExt<ParametersType, DcaDataType, ScalarType>::load_G_II(
    func::function<std::complex<ScalarType>, TotalG4Dmn>& G_II) {
  profiler_type prof(__FUNCTION__, "BseClusterSolverExt", __LINE__);

  if (concurrency.id() == concurrency.first())
    std::cout << "\t" << __FUNCTION__ << " -- over " << KExDmn::dmn_size() << " k exchanges and "
              << WExDmn::dmn_size() << " frequency exchanges.\n\n";

  //std::vector<int> coor_1(G_II.signature(), 0.0);
  // ideally a compile time leaf domain count would be avaialbe from function 
  constexpr int num_indexes_G4 = 10;
  std::array<int, num_indexes_G4> coor_1;
  coor_1.fill(0);
  // Get the first object/channel of the G4 container.
  auto& G4 = data_.get_G4()[0];
  std::array<int, num_indexes_G4> coor_2;
  coor_2.fill(0);
  int ce_size = DCA_matrix_dmn_t::dmn_size();
  std::size_t linind = 0;
  for (int wex_ind = 0; wex_ind < WExDmn::dmn_size(); wex_ind++) {
    for (int kex_ind = 0; kex_ind < KExDmn::dmn_size(); kex_ind++) {
      for (int ce_ind = 0; ce_ind < ce_size; ce_ind++) {
        G4.linind_2_subind(linind++, coor_2);
        // coordinate  0 1 2 3 4 5 6 7 8    9    10 11 12 13 14 15 16   17
        // G4        : b b b b k w k w k_ex w_ex
        // G_II      : b b k w b b k w k_ex w_ex

        coor_1[0] = coor_2[0];
        coor_1[1] = coor_2[1];
        coor_1[2] = coor_2[4];  // k_1
        coor_1[3] = coor_2[5];  // w_1
        coor_1[4] = coor_2[2];
        coor_1[5] = coor_2[3];
        coor_1[6] = coor_2[6];  // k_2
        coor_1[7] = coor_2[7];  // w_2
        coor_1[8] = coor_2[8];  // kex
        coor_1[9] = coor_2[9];  // wex

        // Note: only the first momentum and frequency exchange (9th and 10th index of G4) are analyzed.
        G_II(coor_1) = G4(coor_2);
      }
    }
  }
  std::cout << "final linind " << linind << '\n';
  std::cout << "G_II(0,0,0,0,0,0,0,0,0,0): " << G_II(0,0,0,0,0,0,0,0,0,0) << '\n';
  std::cout << "G_II(0,0,0,0,0,0,0,0,0,1): " << G_II(0,0,0,0,0,0,0,0,0,1) << '\n';
  std::cout << "G_II(0,0,0,0,0,0,0,0,0,2): " << G_II(0,0,0,0,0,0,0,0,0,2) << '\n';
  std::cout << "G4(0,0,0,0,0,0,0,0,0,1): " << G4(0,0,0,0,0,0,0,0,0,1) << '\n';
  std::cout << "G4(0,0,0,0,0,0,0,0,0,2): " << G4(0,0,0,0,0,0,0,0,0,2) << '\n';

}

/** This calculates close to the Chi_0(Q,omega_m)K,K' that we want for S(q,omega)
 */
template <typename ParametersType, typename DcaDataType, typename ScalarType>
void BseClusterSolverExt<ParametersType, DcaDataType, ScalarType>::load_G_II_0(
    func::function<std::complex<ScalarType>, TotalG4Dmn>& G_II_0) {
  profiler_type prof(__FUNCTION__, "BseClusterSolverExt", __LINE__);
  if (concurrency.id() == concurrency.first())
    std::cout << "\t" << __FUNCTION__ << "\n\n";

  func::dmn_variadic<k_DCA, WVertexDmn> k_w_dmn;
  G_II_0 = 0.;

  int W = parameters.get_sp_fermionic_frequencies();  // TODO: Replace using w::dmn_size().
  int W_vertex = WVertexDmn::dmn_size() / 2;
  int coor[2];

  for (int wex_ind = 0; wex_ind < WExDmn::dmn_size(); wex_ind++)
    // Assume all transfers for momentum
    for (int kex_ind = 0; kex_ind < KExDmn::dmn_size(); kex_ind++) {
      for (int i = 0; i < k_w_dmn.get_size(); i++) {
        k_w_dmn.linind_2_subind(i, coor);

        int k = coor[0];
        int w_vertex = coor[1];
        int w = (coor[1] - W_vertex) + W;
        int q = kex_ind;
        int k_plus_q = k_DCA::parameter_type::add(k, q);
        int q_minus_k = k_DCA::parameter_type::subtract(k, q);

        for (int n1 = 0; n1 < b::dmn_size(); n1++) {
          for (int n2 = 0; n2 < b::dmn_size(); n2++) {
            for (int m1 = 0; m1 < b::dmn_size(); m1++) {
              for (int m2 = 0; m2 < b::dmn_size(); m2++) {
                // TODO: allow more than one channel.
                switch (parameters.get_four_point_channels()[0]) {
                  case FourPointType::PARTICLE_HOLE_TRANSVERSE: {
                    G_II_0(n1, n2, k, w_vertex, m1, m2, k, w_vertex, kex_ind, wex_ind) =
                        -data_.G_k_w(n1, e_UP, m2, e_UP, k, w) *
                        data_.G_k_w(n2, e_UP, m1, e_UP, k_plus_q, w + wex_ind);
                    break;
                  }

                  case FourPointType::PARTICLE_HOLE_MAGNETIC: {
                    G_II_0(n1, n2, k, w_vertex, m1, m2, k, w_vertex, 0, wex_ind) =
                        -data_.G_k_w(n1, e_UP, m2, e_UP, k, w) *
                        data_.G_k_w(n2, e_UP, m1, e_UP, k_plus_q, w + wex_ind);
                    break;
                  }
                  case FourPointType::PARTICLE_HOLE_CHARGE: {
                    G_II_0(n1, n2, k, w_vertex, m1, m2, k, w_vertex, 0, wex_ind) =
                        -2. * data_.G_k_w(n1, e_UP, m1, e_UP, k, w) *
                        data_.G_k_w(n2, e_UP, m2, e_UP, k_plus_q, w + wex_ind);
                    break;
                  }

                  case FourPointType::PARTICLE_PARTICLE_UP_DOWN: {
                    double wn = w::get_elements()[w];
                    double w_nu_min_wn = w::get_elements()[wex_ind + (2 * W - 1 - w)];
                    double beta = parameters.get_beta();

                    if (std::fabs((wex_ind * M_PI / beta - wn) - w_nu_min_wn) > 1.e-6)
                      throw std::logic_error(__FUNCTION__);

                    G_II_0(n1, n2, k, w_vertex, m1, m2, k, w_vertex, 0, wex_ind) =
                        data_.G_k_w(n1, e_UP, m1, e_UP, k, w) *
                        data_.G_k_w(n2, e_UP, m2, e_UP, q_minus_k, wex_ind + (2 * W - 1 - w));
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
}

template <typename ParametersType, typename DcaDataType, typename ScalarType>
void BseClusterSolverExt<ParametersType, DcaDataType, ScalarType>::load_G_II_0_function(
    func::function<std::complex<ScalarType>, TotalG4Dmn>& G_II_0) {
  profiler_type prof(__FUNCTION__, "BseClusterSolverExt", __LINE__);
  if (concurrency.id() == concurrency.first())
    std::cout << "\t" << __FUNCTION__ << "\n\n";

  for (int wex_ind = 0; wex_ind < WExDmn::dmn_size(); wex_ind++)
    for (int kex_ind = 0; kex_ind < KExDmn::dmn_size(); kex_ind++)
      for (int w_ind = 0; w_ind < WVertexDmn::dmn_size(); w_ind++)
        for (int K_ind = 0; K_ind < k_DCA::dmn_size(); K_ind++)

          for (int m2 = 0; m2 < b::dmn_size(); m2++)
            for (int n2 = 0; n2 < b::dmn_size(); n2++)

              for (int m1 = 0; m1 < b::dmn_size(); m1++)
                for (int n1 = 0; n1 < b::dmn_size(); n1++)
                  G_II_0_function(n1, m1, n2, m2, K_ind, w_ind, kex_ind, wex_ind) =
                      G_II_0(n1, m1, K_ind, w_ind, n2, m2, K_ind, w_ind, kex_ind, wex_ind);
}

template <typename ParametersType, typename DcaDataType, typename ScalarType>
void BseClusterSolverExt<ParametersType, DcaDataType, ScalarType>::solve_BSE_on_cluster(
    func::function<std::complex<ScalarType>, TotalG4Dmn>& G_II,
    func::function<std::complex<ScalarType>, TotalG4Dmn>& G_II_0) {
  profiler_type prof(__FUNCTION__, "BseClusterSolverExt", __LINE__);
  std::vector<int> subind(3);

  if (concurrency.id() == concurrency.first())
    std::cout << '\t' << __FUNCTION__ << '\n';

  for (int wex_ind = 0; wex_ind < WExDmn::dmn_size(); wex_ind++)
    for (int kex_ind = 0; kex_ind < KExDmn::dmn_size(); kex_ind++) {

      func::function<std::complex<ScalarType>, DCA_matrix_dmn_t> G_II_indi;
      func::function<std::complex<ScalarType>, DCA_matrix_dmn_t> G_II_0_indi;
      subind = {0, kex_ind, wex_ind};
      G_II.slice(0, subind, G_II_indi.values());
      G_II_0.slice(0, subind, G_II_0_indi.values());

      profiler_t prof(__FUNCTION__, __FILE__, __LINE__);

      ScalarType renorm = 1. / (parameters.get_beta() * k_DCA::dmn_size());

      int N = cluster_eigenvector_dmn.get_size();

      G_II_indi *= renorm;
      G_II_0_indi *= renorm;

      // As I understand it each of the below is independent so we should be solving the matrices of the inner bbkwbbkw blocks.
      dca::linalg::Matrix<std::complex<ScalarType>, dca::linalg::CPU> G4_inv_indi(N);
      dca::linalg::Matrix<std::complex<ScalarType>, dca::linalg::CPU> G4_0_inv_indi(N);

      dca::linalg::matrixop::copyArrayToMatrix(N, N, G_II_indi.values(), N, G4_inv_indi);
      dca::linalg::matrixop::inverse(G4_inv_indi);
      dca::linalg::matrixop::copyArrayToMatrix(N, N, G_II_0_indi.values(), N, G4_0_inv_indi);
      dca::linalg::matrixop::inverse(G4_0_inv_indi);
      for (int j = 0; j < N; j++)
        for (int i = 0; i < N; i++)
          Gamma_cluster(i, j, kex_ind, wex_ind) = G4_0_inv_indi(i, j) - G4_inv_indi(i, j);
    }
}
}  // namespace analysis
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_ANALYSIS_BSE_SOLVER_BSE_CLUSTER_SOLVER_EXT_HPP
