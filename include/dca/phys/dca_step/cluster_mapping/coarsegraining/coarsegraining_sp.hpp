// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Giovani Balduzzi(gbalduzz@itp.phys.ethz.ch)
//
// This class performs the coarsegraining of single-particle functions.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_COARSEGRAINING_SP_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_COARSEGRAINING_SP_HPP

#include <complex>
#include <functional>
#include <iostream>
#include <sstream>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/geometry/gaussian_quadrature/gaussian_quadrature_domain.hpp"
#include "dca/math/geometry/tetrahedron_mesh/tetrahedron_mesh.hpp"
#include "dca/phys/dca_algorithms/compute_greens_function.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/coarsegraining_routines.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/interpolation_matrices.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/tetrahedron_integration.hpp"
#include "dca/phys/dca_step/cluster_mapping/coarsegraining/tetrahedron_routines_harmonic_function.hpp"
#include "dca/phys/dca_step/lattice_mapping/interpolation/transform_to_alpha.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_operations.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_domain.hpp"
#include "dca/util/print_time.hpp"
#include "dca/util/plot.hpp"

namespace dca {
namespace phys {
namespace clustermapping {
// dca::phys::clustermapping::

template <typename Parameters>
class CoarsegrainingSp : private coarsegraining_routines<Parameters>,
                         private tetrahedron_integration<Parameters> {
public:
  using Concurrency = typename Parameters::concurrency_type;
  using Threading = typename Parameters::ThreadingType;
  using ThisType = CoarsegrainingSp<Parameters>;

  constexpr static int dimension = Parameters::lattice_dimension;

  using CDA = ClusterDomainAliases<dimension>;
  using KClusterDmn = typename CDA::KClusterDmn;
  using KLatticeDmn = typename CDA::KSpHostDmn;

  using ScalarType = double;
  using Complex = std::complex<ScalarType>;

  using NuDmn =
      func::dmn_variadic<func::dmn_0<domains::electron_band_domain>,
                         func::dmn_0<domains::electron_spin_domain>>;  // orbital-spin index
private:
  using QDmn = func::dmn_0<coarsegraining_domain<KClusterDmn, K>>;

  using WDmn = func::dmn_0<domains::frequency_domain>;

  using ClusterFreqFunction =
      func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn, WDmn>>;
  using LatticeFreqFunction =
      func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KLatticeDmn, WDmn>>;
  using LatticeFunction = func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KLatticeDmn>>;

  using TetDmn = func::dmn_0<coarsegraining_domain<KClusterDmn, TETRAHEDRON_K>>;
  using NuNuDmn = func::dmn_variadic<NuDmn, NuDmn, QDmn>;
  using NuNuTetDmn = func::dmn_variadic<NuDmn, NuDmn, TetDmn>;

public:
  CoarsegrainingSp(Parameters& parameters_ref);

  // Computes the coarsegrained dispersion,
  //     \bar{\varepsilon}_\mathbf{K} = \frac{1}{V} \sum_{\tilde{\mathbf{k}}}
  //                                    \varepsilon_{\mathbf{K}+\tilde{\mathbf{k}}} .
  void computeCoarsegrainedDispersion(
      func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn>>& eps_K_cg) const;

  // Computes the coarse-grained G(K, w) using H_0 and chemichal potential stored in the parameters.
  template <class SigmaType,
            typename = std::enable_if_t<std::is_same<SigmaType, ClusterFreqFunction>::value ||
                                        std::is_same<SigmaType, LatticeFreqFunction>::value>>
  void compute_G_K_w(const SigmaType& S_K_w, ClusterFreqFunction& G_K_w);

  // DCA version of the Sigma coarse-graining. It simply copies S_k_w into S_K_w.
  void compute_S_K_w(const ClusterFreqFunction& S_k_w, ClusterFreqFunction& S_K_w) const;

  // Coarse-grains S_k_w into S_K_w. Used by the  DCA+ algorithm.
  void compute_S_K_w(const LatticeFreqFunction& S_k_w, ClusterFreqFunction& S_K_w);

  // Used exclusively by the analysis application with the DCA+ algorithm.
  template <typename RDmn>
  void compute_phi_r(func::function<ScalarType, RDmn>& phi_r) const;

private:
  template <class SigmaType,
            typename = std::enable_if_t<!std::is_same<SigmaType, LatticeFreqFunction>::value>>
  void updateSigmaInterpolated(const SigmaType& /*Sigma*/) const {}
  void updateSigmaInterpolated(const LatticeFreqFunction& Sigma);

private:
  Parameters& parameters_;
  Concurrency& concurrency_;

  std::vector<func::function<std::complex<ScalarType>, NuNuDmn>> H0_q_;

  // gaussian q-points
  func::function<ScalarType, QDmn> w_q_;
  double w_tot_;

  using SigmaInterpolatedType =
      func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, QDmn, KClusterDmn, WDmn>>;
  std::unique_ptr<SigmaInterpolatedType> Sigma_interpolated_;
  LatticeFreqFunction Sigma_old_;
};

template <typename Parameters>
CoarsegrainingSp<Parameters>::CoarsegrainingSp(Parameters& parameters_ref)
    : coarsegraining_routines<Parameters>(parameters_ref),
      tetrahedron_integration<Parameters>(parameters_ref),

      parameters_(parameters_ref),
      concurrency_(parameters_.get_concurrency()),

      H0_q_(KClusterDmn::dmn_size()),

      w_q_("w_q_"),
      w_tot_(0.) {
  interpolation_matrices<ScalarType, KClusterDmn, QDmn>::initialize(concurrency_);

  // Compute H0(k+q) for each value of k and q.
  for (int k = 0; k < H0_q_.size(); ++k) {
    QDmn::parameter_type::set_elements(k);
    Parameters::model_type::initialize_H_0(parameters_, H0_q_[k]);
  }

  for (int l = 0; l < w_q_.size(); ++l)
    w_tot_ += w_q_(l) = QDmn::parameter_type::get_weights()[l];
}

template <typename Parameters>
void CoarsegrainingSp<Parameters>::computeCoarsegrainedDispersion(
    func::function<Complex, func::dmn_variadic<NuDmn, NuDmn, KClusterDmn>>& eps_K_cg) const {
  eps_K_cg = 0.;

  const func::dmn_variadic<KClusterDmn> K_dmn;
  const std::pair<int, int> concurrency_bounds = concurrency_.get_bounds(K_dmn);

  const int n_threads = parameters_.get_coarsegraining_threads();

  Threading().execute(n_threads, [&](int id, int n_threads) {
    const auto thread_bounds = parallel::util::getBounds(id, n_threads, concurrency_bounds);

    for (int k = thread_bounds.first; k < thread_bounds.second; ++k) {
      const auto& eps_q = H0_q_[k];
      constexpr int n_spin_bands = Parameters::bands * 2;

      for (int q = 0; q < QDmn::dmn_size(); ++q) {
        for (int j = 0; j < n_spin_bands; ++j) {
          for (int i = 0; i < n_spin_bands; ++i) {
            eps_K_cg(i, j, k) += eps_q(i, j, q) * w_q_(q);
          }
        }
      }
    }
  });

  concurrency_.sum(eps_K_cg);

  eps_K_cg /= w_tot_;
}
template <typename Parameters>
template <class SigmaType, typename>
void CoarsegrainingSp<Parameters>::compute_G_K_w(const SigmaType& S_K_w, ClusterFreqFunction& G_K_w) {
  // Computes G_K_w(k,w) = 1/N_q \sum_q 1/(i w + mu - H0(k+q,w) - Sigma(k+q,w)).
  updateSigmaInterpolated(S_K_w);
  G_K_w = 0.;

  func::dmn_variadic<KClusterDmn, WDmn> K_wm_dmn;
  const std::pair<int, int> external_bounds = concurrency_.get_bounds(K_wm_dmn);

  const int n_threads = parameters_.get_coarsegraining_threads();
  Threading().execute(n_threads, [&](int id, int n_threads) {
    const auto bounds = parallel::util::getBounds(id, n_threads, external_bounds);
    linalg::Matrix<Complex, linalg::CPU> G_inv("G_inv", NuDmn::dmn_size());
    linalg::Vector<int, linalg::CPU> ipiv;
    linalg::Vector<Complex, linalg::CPU> work;
    int coor[2];
    func::dmn_variadic<KClusterDmn, WDmn> K_wm_dmn;
    const Complex im(0., 1.);

    for (int l = bounds.first; l < bounds.second; l++) {
      K_wm_dmn.linind_2_subind(l, coor);
      const int k(coor[0]), w(coor[1]);

      const auto w_val = WDmn::get_elements()[w];
      const auto& H0 = H0_q_[k];
      constexpr int n_spin_bands = Parameters::bands * 2;

      for (int q = 0; q < QDmn::dmn_size(); ++q) {
        for (int j = 0; j < n_spin_bands; j++) {
          for (int i = 0; i < n_spin_bands; i++) {
            if (std::is_same<SigmaType, ClusterFreqFunction>::value)
              G_inv(i, j) = -H0(i, j, q) - S_K_w(i, j, k, w);
            else
              G_inv(i, j) = -H0(i, j, q) - (*Sigma_interpolated_)(i, j, q, k, w);
            if (i == j)
              G_inv(i, j) += im * w_val + parameters_.get_chemical_potential();
          }
        }

        linalg::matrixop::inverse(G_inv, ipiv, work);
        for (int j = 0; j < n_spin_bands; ++j)
          for (int i = 0; i < n_spin_bands; ++i)
            G_K_w(i, j, k, w) += G_inv(i, j) * w_q_(q);
      }
    }
  });

  concurrency_.sum(G_K_w);

  G_K_w /= w_tot_;
}

template <typename Parameters>
void CoarsegrainingSp<Parameters>::compute_S_K_w(const ClusterFreqFunction& S_k_w,
                                                 ClusterFreqFunction& S_K_w) const {
  S_K_w = S_k_w;
}

template <typename Parameters>
void CoarsegrainingSp<Parameters>::compute_S_K_w(const LatticeFreqFunction& S_k_w,
                                                 ClusterFreqFunction& S_K_w) {
  S_K_w = 0.;
  updateSigmaInterpolated(S_k_w);

  func::dmn_variadic<KClusterDmn, WDmn> K_wm_dmn;
  const std::pair<int, int> external_bounds = concurrency_.get_bounds(K_wm_dmn);

  Threading().execute(
      parameters_.get_coarsegraining_threads(), [&](const int id, const int n_threads) {
        const auto bounds = parallel::util::getBounds(id, n_threads, external_bounds);

        int coor[2];
        for (int l = bounds.first; l < bounds.second; l++) {
          K_wm_dmn.linind_2_subind(l, coor);
          const int k_ind = coor[0], w_ind = coor[1];

          for (int q_ind = 0; q_ind < QDmn::dmn_size(); q_ind++)
            for (int j = 0; j < NuDmn::dmn_size(); j++)
              for (int i = 0; i < NuDmn::dmn_size(); i++)
                S_K_w(i, j, k_ind, w_ind) +=
                    (*Sigma_interpolated_)(i, j, q_ind, k_ind, w_ind) * w_q_(q_ind);
        }
      });

  concurrency_.sum(S_K_w);

  S_K_w /= w_tot_;
}

template <typename Parameters>
void CoarsegrainingSp<Parameters>::updateSigmaInterpolated(const LatticeFreqFunction& Sigma) {
  if (!Sigma_interpolated_)
    Sigma_interpolated_ = std::make_unique<SigmaInterpolatedType>("Sigma interpolated.");

  if (Sigma_old_ == Sigma)
    return;

  *Sigma_interpolated_ = 0;

  // Compute the interpolation.
  func::dmn_variadic<KClusterDmn, WDmn> K_wm_dmn;
  const std::pair<int, int> external_bounds = concurrency_.get_bounds(K_wm_dmn);
  const int n_threads = parameters_.get_coarsegraining_threads();

  Threading().execute(n_threads, [&](const int id, const int n_threads) {
    const auto bounds = parallel::util::getBounds(id, n_threads, external_bounds);
    LatticeFunction S_k;
    func::function<std::complex<ScalarType>, func::dmn_variadic<NuDmn, NuDmn, QDmn>> S_q;

    int coor[2];
    for (int l = bounds.first; l < bounds.second; l++) {
      K_wm_dmn.linind_2_subind(l, coor);
      const int k(coor[0]), w(coor[1]);
      std::copy_n(&Sigma(0, 0, 0, w), S_k.size(), S_k.values());

      const double alpha = WDmn::get_elements()[w] > 0 ? 1 : -1;
      coarsegraining_routines<Parameters>::wannierInterpolationWithAlphaTransform(k, alpha, S_k, S_q);

      std::copy_n(S_q.values(), S_q.size(), &(*Sigma_interpolated_)(0, 0, 0, k, w));
    }
  });

  concurrency_.sum(*Sigma_interpolated_);

  Sigma_old_ = Sigma;
}

template <typename Parameters>
template <typename RDmn>
void CoarsegrainingSp<Parameters>::compute_phi_r(func::function<ScalarType, RDmn>& phi_r) const {
  using KCluster = typename KClusterDmn::parameter_type;
  math::geometry::tetrahedron_mesh<KCluster> mesh(parameters_.get_k_mesh_recursion());

  using tetrahedron_dmn = func::dmn_0<math::geometry::tetrahedron_mesh<KClusterDmn>>;
  using quadrature_dmn = math::geometry::gaussian_quadrature_domain<tetrahedron_dmn>;
  quadrature_dmn::translate_according_to_period(parameters_.get_coarsegraining_periods(), mesh);

  std::vector<math::geometry::tetrahedron<dimension>>& tetrahedra = mesh.get_tetrahedra();

  phi_r = 0.;

  RDmn r_domain;
  std::pair<int, int> bounds = concurrency_.get_bounds(r_domain);

  std::vector<std::vector<double>> super_basis = RDmn::parameter_type::get_super_basis_vectors();

  for (int l = bounds.first; l < bounds.second; l++) {
    std::vector<double> r_vec = RDmn::get_elements()[l];
    std::vector<std::vector<double>> r_vecs =
        domains::cluster_operations::equivalent_vectors(r_vec, super_basis);
    for (int r_ind = 0; r_ind < r_vecs.size(); r_ind++)
      for (int tet_ind = 0; tet_ind < tetrahedra.size(); tet_ind++)
        phi_r(l) += std::real(tetrahedron_routines_harmonic_function::execute(
                        r_vecs[0], tetrahedra[tet_ind])) /
                    r_vecs.size();
  }

  concurrency_.sum(phi_r);

  double tot_weight = 0;
  for (auto w : QDmn::parameter_type::get_weights())
    tot_weight += w;

  phi_r /= tot_weight;
}

}  // clustermapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_COARSEGRAINING_COARSEGRAINING_SP_HPP
