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
#include <stdexcept>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/linalg.hpp"
#include "dca/math/util/comparison_methods.hpp"
#include "dca/math/util/vector_operations.hpp"
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

  // Number of leading eigenvalues/eigenvectors to store.
  static constexpr int num_evals = 10;
  using LeadingEigDmn = func::dmn_0<func::dmn<num_evals, int>>;

  // Number of cubic harmonics to compare the leading eigenvectors with.
  const static int num_harmonics = 3;
  using CubicHarmonicsDmn = func::dmn_0<func::dmn<num_harmonics, int>>;

  using KExDmn = typename func::dmn_0<domains::MomentumExchangeDomain>;
  using WExDmn = typename func::dmn_0<domains::FrequencyExchangeDomain>;

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
  using crystal_harmonics_expansion = domains::centered_cluster_domain<host_vertex_r_cluster_type>;
  using crystal_harmonics_expansion_dmn_t = func::dmn_0<crystal_harmonics_expansion>;

  using chi_vector_dmn_t = func::dmn_variadic<b, b, crystal_harmonics_expansion_dmn_t>;

  using LatticeEigenvectorDmn = func::dmn_variadic<b, b, k_HOST_VERTEX, WVertexDmn>;

  using crystal_eigenvector_dmn_t =
      func::dmn_variadic<b, b, crystal_harmonics_expansion_dmn_t, WVertexDmn>;
  using CubicHarmonicsEigenvectorDmn = func::dmn_variadic<b, b, CubicHarmonicsDmn, WVertexDmn>;

  using HOST_matrix_dmn_t = func::dmn_variadic<LatticeEigenvectorDmn, LatticeEigenvectorDmn>;

  using SharedDmn = func::dmn_variadic<b, b, WVertexDmn>;
  using SharedMatrixDmn = func::dmn_variadic<SharedDmn, SharedDmn>;
  using GammaLatticeDmn = func::dmn_variadic<b, b, WVertexDmn>;

  using Chi0LatticeDmn = func::dmn_variadic<b_b, b_b, k_HOST_VERTEX, WVertexDmn>;

  using Chi0LatticeMatrixDmn =
      func::dmn_variadic<func::dmn_variadic<b, b, WVertexDmn, b, b, WVertexDmn>, k_HOST_VERTEX, WExDmn>;

  /** A matrix domain over the ChiLatticeDmn, for single site we are only interest in the block
   * diagonal elements where k_HOST, and WVertexDmn are diagonal. */
  using Chi0LatticeMatrixFullDmn =
      func::dmn_variadic<b, b, k_HOST, WVertexDmn, b, b, k_HOST, WVertexDmn>;
  using GLatticeDmn = func::dmn_variadic<b, b, WVertexDmn>;
  /**  results in a sparse matrix that would require sparse linear algebra or extra work to use with dense linalg
   */
  // using Chi0LatticeMatrixDmn = func::dmn_variadic<ChiLatticeDmnb_b, b_b, k_HOST

  using MultExHostLatticeDmn = func::dmn_variadic<HOST_matrix_dmn_t, WExDmn>;

  using ChiDblPrimeDmn = func::dmn_variadic<k_HOST, WExDmn>;

  BseLatticeSolverExt(ParametersType& parameters, DcaDataType& dca_data);

  template <typename Writer>
  void write(Writer& writer);

  void computeChi0Lattice();
  template <typename ClusterEigenDmn, typename KExDmn, typename WExDmn>
  void computeGammaLattice(
      /*const*/ func::function<std::complex<ScalarType>,
                               func::dmn_variadic<ClusterEigenDmn, ClusterEigenDmn, KExDmn, WExDmn>>&
          Gamma_cluster);

  void computeG4Lattice();
  void computeChiDblPrime_q_w();
  void diagonalizeGammaChi0();

  auto& get_leading_eigenvalues() /*const*/ {
    return leading_eigenvalues;
  };
  auto& get_leading_eigenvectors() /*const*/ {
    return leading_eigenvectors;
  };
  auto& get_leading_symmetry_decomposition() /*const*/ {
    return leading_symmetry_decomposition;
  };

private:
  void initialize();

  void diagonalizeGammaChi0Symmetric();
  void diagonalizeGammaChi0Full();
  void diagonalize_folded_Gamma_chi_0();

  template <typename EvElementType>  // Element type of eigenvalues and eigenvectors.
  void recordEigenvaluesAndEigenvectors(const linalg::Vector<EvElementType, linalg::CPU>& L,
                                        const linalg::Matrix<EvElementType, linalg::CPU>& VR);

  ParametersType& parameters;
  concurrency_type& concurrency;

  DcaDataType& dca_data_;

  // only supporting single site
  func::function<std::complex<ScalarType>,
                 func::dmn_variadic<func::dmn_variadic<GammaLatticeDmn, GammaLatticeDmn>, WExDmn>>
      Gamma_lattice;
  func::function<std::complex<ScalarType>, func::dmn_variadic<Chi0LatticeDmn, WExDmn>> chi_0_lattice;
  // Matrix in \vec{k} and \omega_n with the diagonal = chi_0_lattice.
  // func::function<std::complex<ScalarType>, Chi0LatticeMatrixDmn> chi_0_lattice_matrix;

  func::function<std::complex<ScalarType>, LeadingEigDmn> leading_eigenvalues;
  func::function<std::complex<ScalarType>, func::dmn_variadic<LeadingEigDmn, CubicHarmonicsEigenvectorDmn>>
      leading_symmetry_decomposition;
  func::function<std::complex<ScalarType>, func::dmn_variadic<LeadingEigDmn, LatticeEigenvectorDmn>>
      leading_eigenvectors;

  func::function<std::complex<ScalarType>, func::dmn_variadic<chi_vector_dmn_t, chi_vector_dmn_t>> chi_q;

  func::function<std::complex<ScalarType>,
                 func::dmn_variadic<k_HOST_VERTEX, crystal_harmonics_expansion_dmn_t>>
      psi_k;
  func::function<std::complex<ScalarType>,
                 func::dmn_variadic<LatticeEigenvectorDmn, crystal_eigenvector_dmn_t>>
      crystal_harmonics;

  func::function<std::complex<ScalarType>,
                 func::dmn_variadic<func::dmn_variadic<SharedDmn, SharedDmn>, k_HOST, WExDmn>>
      G4_lattice;

  func::function<std::complex<ScalarType>, func::dmn_variadic<k_HOST_VERTEX, CubicHarmonicsDmn>>
      leading_symmetry_functions;

  func::function<std::complex<ScalarType>, ChiDblPrimeDmn> chi_dbl_prime;
};

template <typename ParametersType, typename DcaDataType, typename ScalarType>
BseLatticeSolverExt<ParametersType, DcaDataType, ScalarType>::BseLatticeSolverExt(
    ParametersType& parameters_ref, DcaDataType& dca_data)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      dca_data_(dca_data),

      Gamma_lattice("Gamma_lattice"),
      chi_0_lattice("chi_0_lattice"),
      // chi_0_lattice_matrix("chi_0_lattice-matrix"),

      leading_eigenvalues("leading-eigenvalues"),
      leading_symmetry_decomposition("leading-symmetry-decomposition"),
      leading_eigenvectors("leading-eigenvectors"),

      chi_q("chi(q)"),

      psi_k("psi_k"),
      crystal_harmonics("crystal_harmonics"),
      G4_lattice("G4_lattice"),

      leading_symmetry_functions("leading-symmetry-functions"),
      chi_dbl_prime("chi_dbl_prime") {
  initialize();
}

template <typename ParametersType, typename DcaDataType, typename ScalarType>
template <typename Writer>
void BseLatticeSolverExt<ParametersType, DcaDataType, ScalarType>::write(Writer& writer) {
  writer.execute(Gamma_lattice);
  writer.execute(chi_0_lattice);
  writer.execute(G4_lattice);
  writer.execute(chi_dbl_prime);
}

#define INITDOMAINOUT(x) std::cout << #x << " size " << x::dmn_size() << '\n'

template <typename ParametersType, typename DcaDataType, typename ScalarType>
void BseLatticeSolverExt<ParametersType, DcaDataType, ScalarType>::initialize() {
  {
    INITDOMAINOUT(HOST_matrix_dmn_t);
    INITDOMAINOUT(k_DCA);
    INITDOMAINOUT(k_HOST);
    INITDOMAINOUT(k_HOST_VERTEX);
    INITDOMAINOUT(WVertexDmn);
    INITDOMAINOUT(Chi0LatticeDmn);
    INITDOMAINOUT(Chi0LatticeMatrixDmn);

    // std::cout << "HOST_matrix_dmn_t size : " << HOST_matrix_dmn_t::dmn_size() << '\n';
    // std::cout << "k_HOST_VERTEX size : " << k_HOST_VERTEX::dmn_size() << '\n';

    double r_cut_off = parameters.get_projection_cut_off_radius();

    crystal_harmonics_expansion::initialize();

    std::vector<std::vector<double>> r_vecs(0);

    for (int l = 0; l < crystal_harmonics_expansion::get_size(); l++)
      if (math::util::l2Norm(crystal_harmonics_expansion::get_elements()[l]) < r_cut_off)
        r_vecs.push_back(crystal_harmonics_expansion::get_elements()[l]);

    sort(r_vecs.begin(), r_vecs.end(), math::util::hasSmallerNorm<double>);

    crystal_harmonics_expansion::get_size() = r_vecs.size();
    crystal_harmonics_expansion::get_elements() = r_vecs;

    if (concurrency.id() == concurrency.first()) {
      std::cout << "\n\n\t crystal-vectors : \n";
      for (int l = 0; l < crystal_harmonics_expansion::get_size(); l++) {
        std::cout << "\t" << l << "\t";
        math::util::print(crystal_harmonics_expansion::get_elements()[l]);
        std::cout << std::endl;
      }
    }
  }

  {
    psi_k.reset();
    crystal_harmonics.reset();

    const std::complex<double> I(0, 1);

    for (int l = 0; l < crystal_harmonics_expansion::get_size(); l++) {
      std::vector<double> r_vec = crystal_harmonics_expansion::get_elements()[l];

      for (int k_ind = 0; k_ind < k_HOST_VERTEX::dmn_size(); k_ind++)
        psi_k(k_ind, l) =
            std::exp(I * math::util::innerProduct(r_vec, k_HOST_VERTEX::get_elements()[k_ind])) /
            std::sqrt(double(k_HOST_VERTEX::dmn_size()));

      for (int w_ind = 0; w_ind < WVertexDmn::dmn_size(); w_ind++)
        for (int k_ind = 0; k_ind < k_HOST_VERTEX::dmn_size(); k_ind++)
          for (int m_ind = 0; m_ind < b::dmn_size(); m_ind++)
            for (int n_ind = 0; n_ind < b::dmn_size(); n_ind++)
              crystal_harmonics(n_ind, m_ind, k_ind, w_ind, n_ind, m_ind, l, w_ind) = psi_k(k_ind, l);
    }
  }

  {
    leading_symmetry_functions.reset();
    leading_symmetry_decomposition.reset();

    for (int l = 0; l < CubicHarmonicsDmn::dmn_size(); l++) {
      std::complex<double> norm = 0;

      for (int k_ind = 0; k_ind < k_HOST_VERTEX::dmn_size(); k_ind++) {
        std::complex<double> value = 0;

        double kx = k_HOST_VERTEX::get_elements()[k_ind][0];
        double ky = k_HOST_VERTEX::get_elements()[k_ind][1];

        switch (l) {
          case 0:  // s-wave
            value = 1;
            break;

          case 1:  // p-wave
            value = cos(kx) + cos(ky);
            break;

          case 2:  // d-wave
            value = cos(kx) - cos(ky);
            break;

          default:
            throw std::logic_error(__FUNCTION__);
        }

        norm += conj(value) * value;

        leading_symmetry_functions(k_ind, l) = value;
      }

      for (int k_ind = 0; k_ind < k_HOST_VERTEX::dmn_size(); k_ind++)
        leading_symmetry_functions(k_ind, l) /= std::sqrt(1.e-16 + real(norm));
    }
  }
}

// template <typename ParametersType, typename DcaDataType, typename ScalarType>
// void BseLatticeSolverExt<ParametersType, DcaDataType, ScalarType>::computeChi0LatticeFull() {
//   func::function<std::complex<ScalarType>, GlatticeDmn> g_lattice;
//   auto omega_vertex = WVertexDmn::get_elements();
//   for(int w_ind = 0; w_ind < WVertexDmn::dmn_size(); ++w_ind)
//     for(int k_ind = 0; k_ind < k_HOST::dmn_size(); ++k_ind)
//       for (int m1 = 0; m1 < b::dmn_size(); m1++)
// 	for (int n1 = 0; n1 < b::dmn_size(); n1++)
// 	  g_lattice(n1, m2, k_ind, w_ind) = (omega_vertex(w_ind)+1.0E-16
// }

// TODO: Add finite-size support?
template <typename ParametersType, typename DcaDataType, typename ScalarType>
void BseLatticeSolverExt<ParametersType, DcaDataType, ScalarType>::computeChi0Lattice() {
  profiler_type prof(__FUNCTION__, "BseLatticeSolver", __LINE__);

  std::vector<int> subind(2);
  if (concurrency.id() == concurrency.first())
    std::cout << "\n" << __FUNCTION__ << std::endl;

  clustermapping::coarsegraining_tp<ParametersType, k_HOST_VERTEX> coarsegraining_tp(parameters);
  latticemapping::lattice_mapping_sp<ParametersType, k_DCA, k_HOST> lattice_map_sp(parameters);

  if (parameters.do_dca_plus() || parameters.doPostInterpolation()) {

    dca_data_.Sigma_lattice_interpolated = 0.;
    dca_data_.Sigma_lattice_coarsegrained = 0.;
    dca_data_.Sigma_lattice = 0.;

    if (parameters.hts_approximation()) {
      clustermapping::CoarsegrainingSp<ParametersType> CoarsegrainingSp(parameters);

      DcaDataType dca_data_hts(parameters);
      dca_data_hts.H_HOST = dca_data_.H_HOST;
      dca_data_hts.H_interactions = dca_data_.H_interactions;

      solver::HighTemperatureSeriesExpansionSolver<dca::linalg::CPU, ParametersType, DcaDataType> hts_solver(
          parameters, dca_data_hts);

      lattice_map_sp.execute_with_HTS_approximation(
          dca_data_hts, hts_solver, CoarsegrainingSp, dca_data_.Sigma,
          dca_data_.Sigma_lattice_interpolated, dca_data_.Sigma_lattice_coarsegrained,
          dca_data_.Sigma_lattice);
    }
    else {
      lattice_map_sp.execute(dca_data_.Sigma, dca_data_.Sigma_lattice_interpolated,
                             dca_data_.Sigma_lattice_coarsegrained, dca_data_.Sigma_lattice);
    }
  }

  func::function<std::complex<ScalarType>, Chi0LatticeDmn> chi_0_lattice_indi;
  
  for (int wex_ind = 0; wex_ind < WExDmn::dmn_size(); ++wex_ind) {
    // DCA+/DCA with post-interpolation: Compute \chi_0 from continuous lattice self-energy.
    if (parameters.do_dca_plus() || parameters.doPostInterpolation()) {
      if (parameters.do_dca_plus())
        coarsegraining_tp.execute(dca_data_.H_HOST, dca_data_.Sigma_lattice, chi_0_lattice_indi);
      else  // do_post_interpolation
        coarsegraining_tp.execute(dca_data_.H_HOST, dca_data_.Sigma_lattice_interpolated,
                                  chi_0_lattice_indi);
      // (Standard) DCA: Compute \chi_0 from cluster self-energy.
    }
    else {
      coarsegraining_tp.execute(dca_data_.H_HOST, dca_data_.Sigma, chi_0_lattice_indi);
    }

    subind = {0, wex_ind};
    chi_0_lattice.distribute(0, subind,
                             static_cast<std::complex<ScalarType>*>(chi_0_lattice_indi.values()));
  }
  // Renormalize and set diagonal \chi_0 matrix.
}

// template <typename ParametersType, typename DcaDataType, typename ScalarType>
// template <typename QBrillouin>
// void BseLatticeSolverExt<ParametersType, DcaDataType, ScalarType>::computeG4Lattice() {
//   for(int wex_ind = 0; wex_ind < WExDmn::dmn_size(); ++wex_ind
//     for(int kex_ind = 0; kex_ind < KExDmn::dmn_size(); ++kex_ind)
//       for (int wm_ind = 0; w_ind < WVertexDmn::dmn_size(); wm_ind++)
//         for (int q_ind = 0; K_ind < k_HOST_VERTEX::dmn_size(); q_ind++)
// 	  for (int wn_ind = 0; w_ind < WVertexDmn::dmn_size(); wn_ind++)
// 	for (int wnp_ind = 0; w_ind < WVertexDmn::dmn_size(); wnp_ind++)
//   	for (int m2 = 0; m2 < b::dmn_size(); m2++)
// 	  for (int n2 = 0; n2 < b::dmn_size(); n2++)
// 	    for (int m1 = 0; m1 < b::dmn_size(); m1++)
// 	      for (int n1 = 0; n1 < b::dmn_size(); n1++) {

// 	    }
// }
template <typename ParametersType, typename DcaDataType, typename ScalarType>
template <typename ClusterEigenDmn, typename KExDmn, typename WExDmn>
void BseLatticeSolverExt<ParametersType, DcaDataType, ScalarType>::computeGammaLattice(
    /*const*/ func::function<std::complex<ScalarType>,
                             func::dmn_variadic<ClusterEigenDmn, ClusterEigenDmn, KExDmn, WExDmn>>&
        Gamma_cluster) {
  profiler_type prof(__FUNCTION__, "BseLatticeSolverExt", __LINE__);

  if (concurrency.id() == concurrency.first())
    std::cout << "\n" << __FUNCTION__ << std::endl;

  std::vector<int> subind(4);
  // WExDmn w_dmn;
  // const std::pair<int, int> w_bounds = parallel::util::getBounds(id, nthreads, w_dmn);

  for (int wex_ind = 0; wex_ind < WExDmn::dmn_size(); ++wex_ind) {
    func::function<std::complex<ScalarType>, func::dmn_variadic<ClusterEigenDmn, ClusterEigenDmn>>
        gamma_cluster_indi;
    // Hard coded for single site.
    subind = {0, 0, 0, wex_ind};
    Gamma_cluster.slice(0, 1, subind,
                        static_cast<std::complex<ScalarType>*>(gamma_cluster_indi.values()));
    func::function<std::complex<ScalarType>, func::dmn_variadic<GammaLatticeDmn, GammaLatticeDmn>>
        gamma_lattice_indi;
    for (int w2 = 0; w2 < WVertexDmn::dmn_size(); ++w2)
      for (int m2 = 0; m2 < b::dmn_size(); m2++)
        for (int n2 = 0; n2 < b::dmn_size(); n2++)
          for (int w1 = 0; w1 < WVertexDmn::dmn_size(); ++w1)
            for (int m1 = 0; m1 < b::dmn_size(); m1++)
              for (int n1 = 0; n1 < b::dmn_size(); n1++)
                gamma_lattice_indi(n1, m1, w1, n2, m2, w2) =
                    gamma_cluster_indi(n1, m1, 0, w1, n2, m2, 0, w2);

    std::vector<int> subind_gl{0, wex_ind};
    Gamma_lattice.distribute(0, subind_gl,
                             static_cast<std::complex<ScalarType>*>(gamma_lattice_indi.values()));
  }
  // int n_threads = 1; //WExDmn::dmn_size();
  // Threading parallelization_obj;
  // parallelization_obj.execute(n_threads, makeGammaLatticeOverWExDmn);
}

// template <typename ParametersType, typename DcaDataType, typename ScalarType>
// void BseLatticeSolverExt<ParametersType, DcaDataType, ScalarType>::diagonalizeGammaChi0() {
//   if (parameters.project_onto_crystal_harmonics()) {
//     diagonalize_folded_Gamma_chi_0();
//   }
//   else {
// #ifndef DCA_ANALYSIS_TEST_WITH_FULL_DIAGONALIZATION
//     // Diagonalize the symmetric matrix \sqrt{\chi_0}\Gamma\sqrt{\chi_0}.
//     // The origin in momentum space has always index = 0.
//     // TODO: loop over multiple channels.
//     if (parameters.get_four_point_channels()[0] == FourPointType::PARTICLE_PARTICLE_UP_DOWN &&
//         parameters.get_four_point_momentum_transfer_index() == 0 &&
//         parameters.get_four_point_frequency_transfer() == 0) {
//       diagonalizeGammaChi0Symmetric();
//     }
//     else
// #endif  // DCA_ANALYSIS_TEST_WITH_FULL_DIAGONALIZATION
//       diagonalizeGammaChi0Full();
//   }

//   characterizeLeadingEigenvectors();
//   printOnShell();
// }

template <typename ParametersType, typename DcaDataType, typename ScalarType>
void BseLatticeSolverExt<ParametersType, DcaDataType, ScalarType>::computeG4Lattice() {
  profiler_type prof(__FUNCTION__, "BseLatticeSolverExt", __LINE__);

  if (concurrency.id() == concurrency.first())
    std::cout << "\n" << __FUNCTION__ << std::endl;

  const ScalarType renorm = 1. / (parameters.get_beta() * k_HOST::dmn_size());

  chi_0_lattice *= renorm;

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
                      chi_0_lattice_matrix(n1, m1, w_ind, n2, m2, w_ind, K_ind, wex_ind) =
                          chi_0_lattice(n1, m1, n2, m2, K_ind, w_ind, wex_ind);
                    }
      };

  int n_threads = std::min(8, WExDmn::dmn_size());
  // Threading parallelization_obj;
  Threading().execute(n_threads, makeChi0Matrix, std::ref(chi_0_lattice),
                      std::ref(chi_0_lattice_matrix));

  int N = SharedDmn::dmn_size();
  dca::linalg::Matrix<std::complex<ScalarType>, dca::linalg::CPU> gamma_lat(N);
  dca::linalg::Matrix<std::complex<ScalarType>, dca::linalg::GPU> gamma_lat_GPU(N);

  func::function<std::complex<ScalarType>, SharedMatrixDmn> chi_0_lattice_matrix_indi;
  func::function<std::complex<ScalarType>, SharedMatrixDmn> gamma_lattice_indi;
  std::vector<int> subind(2);

  dca::linalg::Matrix<std::complex<ScalarType>, dca::linalg::CPU> chi0_lat_inv(N);
  dca::linalg::Matrix<std::complex<ScalarType>, dca::linalg::GPU> chi0_lat_inv_GPU(N);

  dca::linalg::Matrix<std::complex<ScalarType>, dca::linalg::CPU> G4_lat_indi(N);
  dca::linalg::Matrix<std::complex<ScalarType>, dca::linalg::GPU> G4_lat_indi_GPU(N);

  for (int wex_ind = 0; wex_ind < WExDmn::dmn_size(); ++wex_ind)
    for (int k_ind = 0; k_ind < k_HOST::dmn_size(); ++k_ind) {
      chi_0_lattice_matrix.slice(
          0, {0, k_ind, wex_ind},
          static_cast<std::complex<ScalarType>*>(chi_0_lattice_matrix_indi.values()));

      dca::linalg::matrixop::copyArrayToMatrix(N, N, chi_0_lattice_matrix_indi.values(), N,
                                               chi0_lat_inv);
      // dca::linalg::util::GpuStream stream;
      chi0_lat_inv_GPU = chi0_lat_inv;
      dca::linalg::matrixop::inverse(chi0_lat_inv_GPU);

      // \todo remove, only necessary to get a look at gamma_lat_inv for debugging.
      chi0_lat_inv = chi0_lat_inv_GPU;

      subind = {0, wex_ind};
      Gamma_lattice.slice(0, subind,
                          static_cast<std::complex<ScalarType>*>(gamma_lattice_indi.values()));

      dca::linalg::matrixop::copyArrayToMatrix(N, N, gamma_lattice_indi.values(), N, gamma_lat);
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
      //   G4_lattice.distribute(0, 1, {0, 0, wex_ind}, G4_lat_inv.ptr());

      for (int j = 0; j < N; j++)
        for (int i = 0; i < N; i++)
          G4_lattice(i + j * N, k_ind, wex_ind) = G4_lat_indi(i, j);
    }
}

template <typename ParametersType, typename DcaDataType, typename ScalarType>
void BseLatticeSolverExt<ParametersType, DcaDataType, ScalarType>::computeChiDblPrime_q_w() {
  profiler_type prof(__FUNCTION__, "BseLatticeSolverExt", __LINE__);

  std::vector<int> g4_lat_subind(3);
  std::size_t shared_matrix_size = SharedDmn::dmn_size() * SharedDmn::dmn_size();
  for (int wex_ind = 0; wex_ind < WExDmn::dmn_size(); ++wex_ind) {
    for (int k_ind = 0; k_ind < k_HOST::dmn_size(); ++k_ind) {
      g4_lat_subind = {0, k_ind, wex_ind};
      func::function<std::complex<ScalarType>, SharedMatrixDmn> g4_lat_indi;
      G4_lattice.slice(0, g4_lat_subind,
                       static_cast<std::complex<ScalarType>*>(g4_lat_indi.values()));
      // auto g4_lat_it = G4_lattice.begin() + G4_lattice.branch_subind_2_linind(g4_lat_subind);
      chi_dbl_prime(k_ind, wex_ind) = std::accumulate(
          g4_lat_indi.begin(), g4_lat_indi.end(), std::complex<ScalarType>{0.0, 0.0}, std::plus<>{});
    }
  }
  ScalarType renorm = 1. / (parameters.get_beta());
  ScalarType renorm_sq = renorm * renorm;
  chi_dbl_prime *= renorm_sq;
}

}  // namespace analysis
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_ANALYSIS_BSE_SOLVER_BSE_LATTICE_SOLVER_HPP
