// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Mi Jiang (jiangmi@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// This class computes the vertex function \Gamma on the lattice and determines the Bethe-Salpeter
// eigenvalues and eigenvectors.

#ifndef DCA_PHYS_DCA_ANALYSIS_BSE_SOLVER_BSE_LATTICE_SOLVER_EXT_HPP
#define DCA_PHYS_DCA_ANALYSIS_BSE_SOLVER_BSE_LATTICE_SOLVER_EXT_HPP

#include <algorithm>
#include <cassert>
#include <complex>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/linalg.hpp"
#include "dca/math/util/comparison_methods.hpp"
#include "dca/math/util/vector_operations.hpp"
#include "dca/phys/dca_analysis/bse_solver/bse_cluster_solver_ext.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/coarsegraining_sp.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/coarsegraining_tp.hpp"
#include "dca/phys/dca_step/cluster_solver/high_temperature_series_expansion/high_temperature_series_expansion_solver.hpp"
#include "dca/phys/dca_step/lattice_mapping/lattice_mapping_sp.hpp"
#include "dca/phys/dca_step/lattice_mapping/lattice_mapping_tp.hpp"
#include "dca/phys/dca_step/symmetrization/diagrammatic_symmetries.hpp"
#include "dca/phys/dca_step/symmetrization/symmetrize.hpp"
#include "dca/phys/domains/cluster/centered_cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"
#include "dca/util/print_time.hpp"

namespace dca {
namespace phys {
namespace analysis {
// dca::phys::analysis::

template <typename ParametersType, typename DcaDataType, typename ScalarType>
class BseLatticeSolverExt {
public:
  using profiler_type = typename ParametersType::profiler_type;
  using concurrency_type = typename ParametersType::concurrency_type;

  using Real = typename util::RealAlias<ScalarType>;

  // Number of leading eigenvalues/eigenvectors to store.
  static constexpr int num_evals = 10;
  using LeadingEigDmn = func::dmn_0<func::dmn<num_evals, int>>;

  using KExDmn = typename func::dmn_0<domains::MomentumExchangeDomain>;
  using WExDmn = typename func::dmn_0<domains::FrequencyExchangeDomain>;

  using WDmn = typename func::dmn_0<domains::frequency_domain>;
  using WVertexDmn = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT>>;
  using b = func::dmn_0<domains::electron_band_domain>;
  using b_b = func::dmn_variadic<b, b>;

  using k_DCA =
      func::dmn_0<domains::cluster_domain<double, ParametersType::lattice_type::DIMENSION, domains::CLUSTER,
                                          domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;
  using k_HOST =
      func::dmn_0<domains::cluster_domain<double, ParametersType::lattice_type::DIMENSION, domains::LATTICE_SP,
                                          domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;
  using k_HOST_VERTEX =
      func::dmn_0<domains::cluster_domain<double, ParametersType::lattice_type::DIMENSION, domains::LATTICE_TP,
                                          domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;

  using host_vertex_r_cluster_type =
      domains::cluster_domain<double, ParametersType::lattice_type::DIMENSION, domains::LATTICE_TP,
                              domains::REAL_SPACE, domains::BRILLOUIN_ZONE>;

  using LatticeEigenvectorDmn = func::dmn_variadic<b, b, k_HOST_VERTEX, WVertexDmn>;

  using HOST_matrix_dmn_t = func::dmn_variadic<LatticeEigenvectorDmn, LatticeEigenvectorDmn>;

  using SharedDmn = func::dmn_variadic<b, b, WVertexDmn>;
  using SharedMatrixDmn = func::dmn_variadic<SharedDmn, SharedDmn>;
  using GammaLatticeDmn = SharedDmn;
  using Chi0LatticeOneQDmn = func::dmn_variadic<SharedDmn, WExDmn>;
  using Chi0LatticeDmn = func::dmn_variadic<Chi0LatticeOneQDmn, k_HOST_VERTEX>;

  using Chi0LatticeOneQ = func::function<std::complex<ScalarType>, Chi0LatticeOneQDmn>;
  using Chi0Lattice = func::function<std::complex<ScalarType>, Chi0LatticeDmn>;

  using OneWExChi0LatticeMatrixDmn =
      func::dmn_variadic<func::dmn_variadic<b, b, WVertexDmn, b, b, WVertexDmn>, k_HOST_VERTEX>;
  using Chi0LatticeMatrixDmn =
      func::dmn_variadic<WExDmn, func::dmn_variadic<b, b, WVertexDmn, b, b, WVertexDmn>, k_HOST_VERTEX>;
  /** A matrix domain over the ChiLatticeDmn, for single site we are only interest in the block
   * diagonal elements where k_HOST, and WVertexDmn are diagonal. */
  using Chi0LatticeMatrixFullDmn =
      func::dmn_variadic<b, b, k_HOST, WVertexDmn, b, b, k_HOST, WVertexDmn>;
  using GLatticeDmn = func::dmn_variadic<b, b, WVertexDmn>;
  /**  results in a sparse matrix that would require sparse linear algebra or extra work to use with dense linalg
   */

  using G4LatticeOneQDmn = func::dmn_variadic<SharedMatrixDmn, WExDmn>;
  using G4LatticeDmn = func::dmn_variadic<G4LatticeOneQDmn, k_HOST>;

  using G4LatticeOneQ = func::function<std::complex<ScalarType>, G4LatticeOneQDmn>;
  using G4Lattice = func::function<std::complex<ScalarType>, G4LatticeDmn>;

  using MultExHostLatticeDmn = func::dmn_variadic<HOST_matrix_dmn_t, WExDmn>;

  using ChiDblPrimeDmn = func::dmn_variadic<k_HOST, WExDmn>;

  BseLatticeSolverExt(ParametersType& parameters, DcaDataType& dca_data);

  template <typename Writer>
  void write(Writer& writer);

  // void computeChi0Lattice();
  void computeChi0LatticeOverHost();
  /** Compute Chi0LatticeSingleSite as in the python miniapp */
  void computeChi0LatticeSingleSite(const std::vector<double>& chi_q, Chi0LatticeOneQ& cloq);

  template <typename DCA_matrix_dmn_t, typename KExDmn, typename WExDmn>
  void computeGammaLattice(
      /*const*/ func::function<std::complex<ScalarType>,
                               func::dmn_variadic<DCA_matrix_dmn_t, KExDmn, WExDmn>>& Gamma_cluster);

  void computeG4Lattice();
  void computeG4LatticeOverHost();
  void computeG4LatticeSingleSite(G4LatticeOneQ& g4l, Chi0LatticeOneQ& chi_0_lattice_one_q);
  void computeChiDblPrime_q_w();
  void diagonalizeGammaChi0();

private:
  void initialize();

  void diagonalizeGammaChi0Symmetric();
  void diagonalizeGammaChi0Full();
  void diagonalize_folded_Gamma_chi_0();

  template <typename EvElementType>  // Element type of eigenvalues and eigenvectors.
  void recordEigenvaluesAndEigenvectors(const linalg::Vector<EvElementType, linalg::CPU>& L,
                                        const linalg::Matrix<EvElementType, linalg::CPU>& VR);

  auto makeDispersionFunc(std::vector<double> q) {
    assert(q.size() == 2);
    double model_t = parameters.get_t();
    double t_prime = parameters.get_t_prime();
    return [model_t, t_prime, q](auto& k_elem) -> double {
      return -2. * model_t * (std::cos(k_elem[0] + q[0]) + std::cos(k_elem[1] + q[1])) -
             4.0 * t_prime * std::cos(k_elem[0] + q[0]) * std::cos(k_elem[1] + q[1]);
    };
  }

  ParametersType& parameters;
  concurrency_type& concurrency;

  DcaDataType& dca_data_;

  // only supporting single site
  func::function<std::complex<ScalarType>, func::dmn_variadic<SharedMatrixDmn, WExDmn>> gamma_lattice_;

  Chi0Lattice chi_0_lattice;

  typename parallel::thread_traits::mutex_type g4l_mutex_;

  G4Lattice G4_lattice;

  func::function<std::complex<ScalarType>, ChiDblPrimeDmn> chi_dbl_prime;

  /// for optimization of fine k grid for chi_0_lattice calculation
  static constexpr int n_k_{100};
  static constexpr int n_k_grid_{n_k_ * n_k_};

  std::vector<double> ek_;
  std::vector<std::array<double, 2>> k_grid_fine_;
};

template <typename ParametersType, typename DcaDataType, typename ScalarType>
BseLatticeSolverExt<ParametersType, DcaDataType, ScalarType>::BseLatticeSolverExt(
    ParametersType& parameters_ref, DcaDataType& dca_data)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),
      dca_data_(dca_data),
      gamma_lattice_("Gamma_lattice"),
      chi_0_lattice("chi_0_lattice"),
      G4_lattice("G4_lattice"),
      chi_dbl_prime("chi_dbl_prime"),
      ek_(n_k_grid_),
      k_grid_fine_(n_k_grid_) {
  initialize();
}

template <typename ParametersType, typename DcaDataType, typename ScalarType>
template <typename Writer>
void BseLatticeSolverExt<ParametersType, DcaDataType, ScalarType>::write(Writer& writer) {
  writer.execute(gamma_lattice_);
  writer.execute(chi_0_lattice);
  writer.execute(G4_lattice);
  writer.execute(chi_dbl_prime);
}

#define INITDOMAINOUT(x) std::cout << #x << " size " << x::dmn_size() << '\n'

template <typename ParametersType, typename DcaDataType, typename ScalarType>
void BseLatticeSolverExt<ParametersType, DcaDataType, ScalarType>::initialize() {
  {
    INITDOMAINOUT(k_DCA);
    INITDOMAINOUT(k_HOST);
    INITDOMAINOUT(k_HOST_VERTEX);
    INITDOMAINOUT(WVertexDmn);
    INITDOMAINOUT(Chi0LatticeDmn);
    INITDOMAINOUT(Chi0LatticeOneQDmn);

    // std::cout << "HOST_matrix_dmn_t size : " << HOST_matrix_dmn_t::dmn_size() << '\n';
    // std::cout << "k_HOST_VERTEX size : " << k_HOST_VERTEX::dmn_size() << '\n';
  }

  // initialize the fine k_grid for chi_0_lattice calculation
  std::vector<double> k_elements(n_k_);
  for (int ik = 0; ik < n_k_; ++ik)
    k_elements[ik] = -M_PI + ((2 * M_PI) / n_k_) * ik;

  auto iter_elements = k_grid_fine_.begin();
  for (int iky = 0; iky < n_k_; ++iky)
    for (int ikx = 0; ikx < n_k_; ++ikx) {
      (*iter_elements) = {k_elements[ikx], k_elements[iky]};
      ++iter_elements;
    }

  std::transform(k_grid_fine_.begin(), k_grid_fine_.end(), ek_.begin(),
                 makeDispersionFunc({0.0, 0.0}));
}

template <typename ParametersType, typename DcaDataType, typename ScalarType>
void BseLatticeSolverExt<ParametersType, DcaDataType, ScalarType>::computeChi0LatticeOverHost() {
  std::vector<std::vector<double>> k_elements = k_HOST::get_elements();
  std::cout << "K_ELEMENTS!\n";
  std::cout << vectorToString(k_elements) << '\n';

  int n_threads = 8;

  std::vector<Chi0LatticeOneQ> chi0_qs(n_threads);

  auto computeSingleSiteChi0 = [&](auto& k_elems, auto& chi0_q) {
    computeChi0LatticeSingleSite(k_elems, chi0_q);
  };

  auto makeChi0LatticeOneQ = [&computeSingleSiteChi0, &k_elements](
                                 int id, int nr_threads, std::vector<Chi0LatticeOneQ>& chi0_qs,
                                 Chi0Lattice& chi_0_lattice) {
    k_HOST k_dmn;
    std::vector<int> subind(2);
    const std::pair<int, int> k_bounds = parallel::util::getBounds(id, nr_threads, k_dmn);
    for (int k_ind = k_bounds.first; k_ind < k_bounds.second; ++k_ind) {
      computeSingleSiteChi0(k_elements[k_ind], chi0_qs[id]);
      subind = {0, k_ind};
      // could cause cache line contention.
      chi_0_lattice.distribute(0, subind,
                               static_cast<std::complex<ScalarType>*>(chi0_qs[id].values()));
      std::cout << '.';
    }
  };

  std::cout << "Calculating Chi0Lattice on " << n_threads << " threads:";

  Threading().execute(n_threads, makeChi0LatticeOneQ, std::ref(chi0_qs), std::ref(chi_0_lattice));
  std::cout << '\n';
}

template <typename ParametersType, typename DcaDataType, typename ScalarType>
void BseLatticeSolverExt<ParametersType, DcaDataType, ScalarType>::computeChi0LatticeSingleSite(
    const std::vector<double>& chi_q, Chi0LatticeOneQ& cloq) {
  profiler_type prof(__FUNCTION__, "BseLatticeSolver", __LINE__);

  std::vector<double> ekpq(n_k_grid_);

  // std::cout << "WVertexDmn::dmn_size()" << WVertexDmn::dmn_size() << '\n';
  int n_w_G4 = WVertexDmn::dmn_size();
  int mid_index_w_G40 = n_w_G4 / 2;

  int n_w_G = WDmn::dmn_size();
  int mid_index_w_G0 = n_w_G / 2;

  // these are a good candidate for memoization
  std::transform(k_grid_fine_.begin(), k_grid_fine_.end(), ekpq.begin(), makeDispersionFunc(chi_q));

  auto& w_set = WDmn::get_elements();
  auto& wn_set = WVertexDmn::get_elements();
  double inv_beta = 1 / parameters.get_beta();
#ifndef NDEBUG
  // auto& wex_set = WExDmn::get_elements();
  // std::cout << "WVertexDmn::elements: " << vectorToString(wn_set) << '\n';
  // std::cout << "WExDmn::elements: " << vectorToString(wex_set) << '\nf';
  // std::vector<double> wex_translated(wex_set.size(), 0.0);
  // std::transform(wex_set.begin(), wex_set.end(), wex_translated.begin(),
  //                [inv_beta](int wex_ind) { return 2 * wex_ind * M_PI * inv_beta; });
  // std::cout << "WEx::elements translated: " << vectorToString(wex_translated) << '\n';
#endif

  double mu = parameters.get_chemical_potential();
  std::complex<double> inv_n_k_grid{1.0 / static_cast<double>(n_k_grid_), 0.0};

  std::vector<std::complex<double>> cc(n_k_grid_);
  for (int wex_ind = 0; wex_ind < WExDmn::dmn_size(); ++wex_ind) {
    for (int iwn = 0; iwn < n_w_G4; ++iwn) {
      // wn is an actual frequency from WVertexDmn
      int i_wG = iwn - mid_index_w_G40 + mid_index_w_G0;
      int iwPlusiwm = std::min(std::max(i_wG + wex_ind, 0), n_w_G - 1);  // iwn + iwm
      double wex = 2 * wex_ind * M_PI * inv_beta;
      assert((wex + wn_set[iwn]) - w_set[iwPlusiwm] < 0.0001);
      for (int ik = 0; ik < n_k_grid_; ++ik) {
        auto& sigma = dca_data_.Sigma;
        std::complex<double> c1{0, wn_set[iwn]};
        std::complex<double> c2{0, wex + wn_set[iwn]};
        // assert(ek[ik] != ekpq[ik]);
        c1 = 1.0 / (c1 + (mu - ek_[ik] - sigma(0, 0, 0, i_wG)));
        c2 = 1.0 / (c2 + (mu - ekpq[ik] - sigma(0, 0, 0, iwPlusiwm)));
        cc[ik] = -c1 * c2;
      }
      std::complex<double> chi0_elem =
          std::accumulate(cc.begin(), cc.end(), std::complex<double>{0, 0});
      cloq(0, 0, iwn, wex_ind) = chi0_elem * inv_n_k_grid;
    }
  }
}

template <typename ParametersType, typename DcaDataType, typename ScalarType>
template <typename ClusterEigenMatrixDmn, typename KExDmn, typename WExDmn>
void BseLatticeSolverExt<ParametersType, DcaDataType, ScalarType>::computeGammaLattice(
    /*const*/ func::function<std::complex<ScalarType>,
                             func::dmn_variadic<ClusterEigenMatrixDmn, KExDmn, WExDmn>>& Gamma_cluster) {
  profiler_type prof(__FUNCTION__, "BseLatticeSolverExt", __LINE__);

  if (concurrency.id() == concurrency.first())
    std::cout << "\n" << __FUNCTION__ << std::endl;

  std::vector<int> subind(4);
  // WExDmn w_dmn;
  // const std::pair<int, int> w_bounds = parallel::util::getBounds(id, nthreads, w_dmn);

  double inv_beta = 1 / parameters.get_beta();

  for (int wex_ind = 0; wex_ind < WExDmn::dmn_size(); ++wex_ind) {
    func::function<std::complex<ScalarType>, ClusterEigenMatrixDmn> gamma_cluster_indi;
    // Hard coded for single site.
    subind = {0, 0, wex_ind};
    Gamma_cluster.slice(0, subind,
                        static_cast<std::complex<ScalarType>*>(gamma_cluster_indi.values()));
    func::function<std::complex<ScalarType>, SharedMatrixDmn> gamma_lattice_indi;
    for (int w2 = 0; w2 < WVertexDmn::dmn_size(); ++w2)
      for (int m2 = 0; m2 < b::dmn_size(); m2++)
        for (int n2 = 0; n2 < b::dmn_size(); n2++)
          for (int w1 = 0; w1 < WVertexDmn::dmn_size(); ++w1)
            for (int m1 = 0; m1 < b::dmn_size(); m1++)
              for (int n1 = 0; n1 < b::dmn_size(); n1++)
                gamma_lattice_indi(n1, m1, w1, n2, m2, w2) =
                    gamma_cluster_indi(n1, m1, 0, w1, n2, m2, 0, w2);

    std::vector<int> subind_gl{0, wex_ind};
    gamma_lattice_.distribute(0, subind_gl,
                              static_cast<std::complex<ScalarType>*>(gamma_lattice_indi.values()));
  }
}

template <typename ParametersType, typename DcaDataType, typename ScalarType>
void BseLatticeSolverExt<ParametersType, DcaDataType, ScalarType>::computeG4Lattice() {
  profiler_type prof(__FUNCTION__, "BseLatticeSolverExt", __LINE__);

  if (concurrency.id() == concurrency.first())
    std::cout << "\n" << __FUNCTION__ << std::endl;

  func::function<std::complex<ScalarType>, Chi0LatticeMatrixDmn> chi_0_lattice_matrix;

  auto makeChi0Matrix =
      [](int id, int nr_threads,
         func::function<std::complex<ScalarType>, func::dmn_variadic<Chi0LatticeDmn, WExDmn>>& chi_0_lattice,
         func::function<std::complex<ScalarType>, Chi0LatticeMatrixDmn>& chi_0_lattice_matrix) {
        WExDmn wex_dmn;
        const std::pair<int, int> wex_bounds = parallel::util::getBounds(id, nr_threads, wex_dmn);
        for (int wex_ind = wex_bounds.first; wex_ind < wex_bounds.second; wex_ind++)
          for (int w_ind = 0; w_ind < WVertexDmn::dmn_size(); w_ind++)
            for (int K_ind = 0; K_ind < k_HOST::dmn_size(); K_ind++)
              for (int m2 = 0; m2 < b::dmn_size(); m2++)
                for (int n2 = 0; n2 < b::dmn_size(); n2++)
                  for (int m1 = 0; m1 < b::dmn_size(); m1++)
                    for (int n1 = 0; n1 < b::dmn_size(); n1++) {
                      // This is block diagonal for single site only
                      chi_0_lattice_matrix(wex_ind, n1, m1, w_ind, n2, m2, w_ind, K_ind, wex_ind) =
                          chi_0_lattice(n1, m1, n2, m2, K_ind, w_ind);
                    }
      };

  int n_threads = 1;  // std::min(8, WExDmn::dmn_size());
  // Threading parallelization_obj;
  Threading().execute(n_threads, makeChi0Matrix, std::ref(chi_0_lattice),
                      std::ref(chi_0_lattice_matrix));

  int N = SharedDmn::dmn_size();
  dca::linalg::Matrix<std::complex<ScalarType>, dca::linalg::CPU> gamma_lat(N);

  func::function<std::complex<ScalarType>, SharedMatrixDmn> chi_0_lattice_matrix_indi;
  func::function<std::complex<ScalarType>, SharedMatrixDmn> gamma_lattice_indi;
  std::vector<int> subind(2);

  dca::linalg::Matrix<std::complex<ScalarType>, dca::linalg::CPU> chi0_lat_inv(N);
  dca::linalg::Matrix<std::complex<ScalarType>, dca::linalg::CPU> G4_lat_indi(N);

#ifdef DCA_HAVE_GPU
  dca::linalg::Matrix<std::complex<ScalarType>, dca::linalg::GPU> gamma_lat_GPU(N);
  dca::linalg::Matrix<std::complex<ScalarType>, dca::linalg::GPU> chi0_lat_inv_GPU(N);
  dca::linalg::Matrix<std::complex<ScalarType>, dca::linalg::GPU> G4_lat_indi_GPU(N);
#endif

  for (int wex_ind = 0; wex_ind < WExDmn::dmn_size(); ++wex_ind)
    for (int k_ind = 0; k_ind < k_HOST::dmn_size(); ++k_ind) {
      chi_0_lattice_matrix.slice(
          0, {0, k_ind, wex_ind},
          static_cast<std::complex<ScalarType>*>(chi_0_lattice_matrix_indi.values()));

      dca::linalg::matrixop::copyArrayToMatrix(N, N, chi_0_lattice_matrix_indi.values(), N,
                                               chi0_lat_inv);
      // dca::linalg::util::GpuStream stream;
#ifdef DCA_HAVE_GPU
      chi0_lat_inv_GPU = chi0_lat_inv;
      dca::linalg::matrixop::inverse(chi0_lat_inv_GPU);

      // \todo remove, only necessary to get a look at gamma_lat_inv for debugging.
      chi0_lat_inv = chi0_lat_inv_GPU;
#else
      dca::linalg::matrixop::inverse(chi0_lat_inv);
#endif
      subind = {0, wex_ind};
      gamma_lattice_.slice(0, subind,
                           static_cast<std::complex<ScalarType>*>(gamma_lattice_indi.values()));

      dca::linalg::matrixop::copyArrayToMatrix(N, N, gamma_lattice_indi.values(), N, gamma_lat);

#ifdef DCA_HAVE_GPU
      gamma_lat_GPU = gamma_lat;
      G4_lat_indi_GPU = chi0_lat_inv_GPU;

      // for (int j = 0; j < N; j++)
      //   for (int i = 0; i < N; i++)
      //     G4_lat_inv(i, j) = chi0_lat_inv(i, j) - gamma_lat_inv(i, j);

      // chi_0_lat_inv_GPU ends up with G4_lat_inv
      dca::linalg::blas::UseDevice<linalg::GPU>::axpy(N * N, std::complex<ScalarType>(-1.0, 0.0),
                                                      gamma_lat_GPU.ptr(), 1, G4_lat_indi_GPU.ptr(),
                                                      1, 0, 0);
      // This is the inv of G4_lat_inv
      dca::linalg::matrixop::inverse(G4_lat_indi_GPU);
      G4_lat_indi = G4_lat_indi_GPU;
#else
      G4_lat_indi = chi0_lat_inv;
      dca::linalg::blas::UseDevice<linalg::CPU>::axpy(N * N, std::complex<ScalarType>(-1.0, 0.0),
                                                      gamma_lat.ptr(), 1, G4_lat_indi.ptr(), 1, 0, 0);
      dca::linalg::matrixop::inverse(G4_lat_indi);
#endif
      for (int j = 0; j < N; j++)
        for (int i = 0; i < N; i++)
          G4_lattice(i + j * N, k_ind, wex_ind) = G4_lat_indi(i, j);
    }
}

template <typename ParametersType, typename DcaDataType, typename ScalarType>
void BseLatticeSolverExt<ParametersType, DcaDataType, ScalarType>::computeG4LatticeOverHost() {
  profiler_type prof(__FUNCTION__, "BseLatticeSolverExt", __LINE__);
  if (concurrency.id() == concurrency.first())
    std::cout << "\n" << __FUNCTION__ << std::endl;

  int num_wex = WExDmn::dmn_size();
  double inv_beta = 1 / parameters.get_beta();

  int n_threads = 8;

  std::vector<Chi0LatticeOneQ> chi0_qs(n_threads);

  std::vector<G4LatticeOneQ> g4ls(n_threads);

  auto computeG4OneQ = [&](G4LatticeOneQ& g4l, Chi0LatticeOneQ& chi0_q) {
    computeG4LatticeSingleSite(g4l, chi0_q);
  };

  auto writeG4OneQ = [&](G4LatticeOneQ& g4l, int k_ind) {
    parallel::thread_traits::scoped_lock lock(g4l_mutex_);
    std::vector<int> subind_g4{0, k_ind};
    std::size_t linind = G4_lattice.branch_subind_2_linind(subind_g4);
    std::copy(g4l.begin(), g4l.end(), G4_lattice.begin() + linind);
    // G4_lattice.distribute(0, subind_g4, static_cast<std::complex<ScalarType>*>(g4l.values()));
  };

  auto makeG4LatticeOneQ = [&computeG4OneQ, &writeG4OneQ](
                               int id, int nr_threads, std::vector<Chi0LatticeOneQ>& chi0_qs,
                               std::vector<G4LatticeOneQ>& g4ls, Chi0Lattice& chi_0_lattice) {
    k_HOST k_dmn;
    const std::pair<int, int> k_bounds = parallel::util::getBounds(id, nr_threads, k_dmn);
    for (int k_ind = k_bounds.first; k_ind < k_bounds.second; ++k_ind) {
      std::vector<int> subind = {0, k_ind};
      chi_0_lattice.slice(0, subind, static_cast<std::complex<ScalarType>*>(chi0_qs[id].values()));
      computeG4OneQ(g4ls[id], chi0_qs[id]);
      writeG4OneQ(g4ls[id], k_ind);
      // mutex bulk assign G4_lattice
    }
  };

  Threading().execute(n_threads, makeG4LatticeOneQ, std::ref(chi0_qs), std::ref(g4ls),
                      std::ref(chi_0_lattice));
}

template <typename ParametersType, typename DcaDataType, typename ScalarType>
void BseLatticeSolverExt<ParametersType, DcaDataType, ScalarType>::computeG4LatticeSingleSite(
    G4LatticeOneQ& g4l, Chi0LatticeOneQ& chi_0_lattice_one_q) {
  profiler_type prof(__FUNCTION__, "BseLatticeSolverExt", __LINE__);

  int num_wex = WExDmn::dmn_size();
  double beta = parameters.get_beta();
  double inv_beta = 1 / parameters.get_beta();

  for (int iwex = 0; iwex < num_wex; ++iwex) {
    dca::linalg::Vector<std::complex<ScalarType>, dca::linalg::CPU> diag_wn(
        chi_0_lattice_one_q.get_domain().get_branch_size(0));
    chi_0_lattice_one_q.slice(0, {0, iwex}, diag_wn.data());

    // Matrix may have aligned rows/cols
    auto chi0_diag = dca::linalg::makeDiagonalMatrixInv(diag_wn);
    int N = diag_wn.size();
    dca::linalg::Matrix<std::complex<ScalarType>, dca::linalg::CPU> g2l_final("g2l", N);
    dca::linalg::Matrix<std::complex<ScalarType>, dca::linalg::CPU> g2l("g2l", N);
    dca::linalg::Matrix<std::complex<ScalarType>, dca::linalg::CPU> gamma_mat("gamma", N);
    func::function<std::complex<ScalarType>, SharedMatrixDmn> gamma_lattice_indi;
    std::size_t linind = gamma_lattice_.branch_subind_2_linind({0, iwex});
    // std::copy(gamma_lattice_.begin() + linind, gamma_lattice_.begin() + linind +
    // gamma_lattice_indi.size(), gamma_lattice_indi.begin()) gamma_lattice_.slice(0, {0, iwex},
    //                      static_cast<std::complex<ScalarType>*>(gamma_lattice_indi.values()));
    // gamma_lattice_indi *= beta;
    dca::linalg::matrixop::copyArrayToMatrix(N, N, gamma_lattice_.values() + linind, N, gamma_mat);
    g2l = chi0_diag;

    assert(gamma_mat.leadingDimension() == g2l.leadingDimension());
    for (int j = 0; j < N; ++j)
      for (int i = 0; i < N; ++i)
        g2l_final(i, j) = g2l(i, j) - gamma_mat(i, j) * inv_beta;

    // for (int ir = 0; ir < gamma_mat.nrRows(); ++ir)
    //   dca::linalg::blas::UseDevice<linalg::CPU>::axpy(
    //       N, std::complex<ScalarType>(-1.0, 0.0), gamma_mat.ptr(ir, 0), N, g2l.ptr(ir, 0), N, 0, 0);

    dca::linalg::matrixop::inverse(g2l_final);

    // we can't copy this back in a faster fashion because it quite likely g2l is not contiguous.
    for (int j = 0; j < N; ++j)
      for (int i = 0; i < N; ++i)
        g4l(0, 0, i, 0, 0, j, iwex) = g2l_final(i, j);
  }
}

template <typename ParametersType, typename DcaDataType, typename ScalarType>
void BseLatticeSolverExt<ParametersType, DcaDataType, ScalarType>::computeChiDblPrime_q_w() {
  profiler_type prof(__FUNCTION__, "BseLatticeSolverExt", __LINE__);

  std::vector<int> g4_lat_subind(8);
  std::size_t shared_matrix_size = SharedDmn::dmn_size() * SharedDmn::dmn_size();
  // func::function<std::complex<ScalarType>, SharedMatrixDmn> g4_lat_indi;
  for (int wex_ind = 0; wex_ind < WExDmn::dmn_size(); ++wex_ind) {
    for (int k_ind = 0; k_ind < k_HOST::dmn_size(); ++k_ind) {
      g4_lat_subind = {0, 0, 0, 0, 0, 0, wex_ind, k_ind};
      std::size_t linind = G4_lattice.subind_2_linind(g4_lat_subind);
      auto g4_lat_begin = G4_lattice.begin() + linind;
      auto g4_lat_end = G4_lattice.begin() + linind + shared_matrix_size;
      // std::copy(G4_lattice.begin()+linind, g4_lat_end, g4_lat_indi.begin());
      //  G4_lattice.slice(5, g4_lat_subind,
      //                   static_cast<std::complex<ScalarType>*>(g4_lat_indi.values()));
      //  auto g4_lat_it = G4_lattice.begin() + G4_lattice.branch_subind_2_linind(g4_lat_subind);
      chi_dbl_prime(k_ind, wex_ind) = std::accumulate(
          g4_lat_begin, g4_lat_end, std::complex<ScalarType>{0.0, 0.0}, std::plus<>{});
    }
  }
  ScalarType renorm = 1. / (parameters.get_beta());
  chi_dbl_prime *= renorm;
}
}  // namespace analysis
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_ANALYSIS_BSE_SOLVER_BSE_LATTICE_SOLVER_HPP
