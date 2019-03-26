// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class computes the lattice self-energy from the dual self-energy.

#ifndef DCA_PHYS_DUAL_FERMION_DUAL_TO_LATTICE_SELF_ENERGY_HPP
#define DCA_PHYS_DUAL_FERMION_DUAL_TO_LATTICE_SELF_ENERGY_HPP

#include <cassert>
#include <complex>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/math/function_transform/special_transforms/space_transform_2D.hpp"
#include "dca/math/util/vector_operations.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/phys/domains/cluster/cluster_operations.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"

namespace dca {
namespace phys {
namespace df {
// dca::phys::df::

template <typename Scalar, typename Concurrency, int dimension>
class DualToLatticeSelfEnergy {
public:
  using Complex = std::complex<Scalar>;

  using BandDmn = func::dmn_0<phys::domains::electron_band_domain>;
  using SpinDmn = func::dmn_0<phys::domains::electron_spin_domain>;

  using SpFreqDmn = func::dmn_0<phys::domains::frequency_domain>;

  using RClusterDmn = typename phys::ClusterDomainAliases<dimension>::RClusterDmn;
  using KClusterDmn = typename phys::ClusterDomainAliases<dimension>::KClusterDmn;
  using KLatticeDmn = typename phys::ClusterDomainAliases<dimension>::KSpHostDmn;
  using KSuperlatticeDmn = typename phys::ClusterDomainAliases<dimension>::KSpSuperlatticeDmn;

  using DualGFDmn = func::dmn_variadic<KClusterDmn, KClusterDmn, KSuperlatticeDmn, SpFreqDmn>;
  using DualGF = func::function<Complex, DualGFDmn>;

  using SpClusterGFDmn =
      func::dmn_variadic<BandDmn, SpinDmn, BandDmn, SpinDmn, KClusterDmn, SpFreqDmn>;
  using SpClusterGF = func::function<Complex, SpClusterGFDmn>;

  using SpLatticeGFDmn =
      func::dmn_variadic<BandDmn, SpinDmn, BandDmn, SpinDmn, KLatticeDmn, SpFreqDmn>;
  using SpLatticeGF = func::function<Complex, SpLatticeGFDmn>;

  DualToLatticeSelfEnergy(const Concurrency& concurrency, const SpClusterGF& G_cluster,
                          const SpClusterGF& Sigma_cluster, const DualGF& Sigma_dual,
                          SpLatticeGF& Sigma_lattice)
      : concurrency_(concurrency),
        G_cluster_(G_cluster),
        Sigma_cluster_(Sigma_cluster),
        Sigma_dual_(Sigma_dual),
        Sigma_lattice_(Sigma_lattice) {
    // TODO: Multi-orbital support.
    assert(BandDmn::dmn_size() == 1);
  }

  // Computes the non-diagonal (in K) lattice self-energy from the dual self-energy,
  //     Sigma(K, K', k_tilde, w) = Sigma_c(K, w) delta_{K, K'} + Sigma_bar(K, K', k_tilde, w),
  // with (matrix notation in K, K')
  //     Sigma_bar = (1 + Sigma_dual * G_c)^{-1} * Sigma_dual.
  void computeNonDiagonalLatticeSelfEnergy();

  // Computes the diagonal lattice self-energy from the non-diagonal (in K) lattice self-energy,
  //     Sigma(k = K + k_tilde) = 1/Nc sum_{I, J} Sigma(I, J, k_tilde) e^{i (I - J) * (K + k_tilde)},
  // with
  //    Sigma(I, J, k_tilde) = 1/Nc sum_{K, K'} Sigma(K, K', k_tilde) e^{-i (I * K - J * K')}.
  void computeDiagonalLatticeSelfEnergy() {}

  const DualGF& getNonDiagonalLatticeSelfEnergy() const {
    return Sigma_lattice_nondiag_;
  }

  static void makeTranslationalInvariant(
      const func::function<Complex, func::dmn_variadic<RClusterDmn, RClusterDmn>>& f_in,
      func::function<Complex, RClusterDmn>& f_out);

  static void updateLatticeSelfEnergy(const func::function<Complex, KClusterDmn>& Sigma_lattice_K,
                                      const int k_tilde, const int w, SpLatticeGF& Sigma_lattice);

private:
  const Concurrency& concurrency_;

  const SpClusterGF& G_cluster_;
  const SpClusterGF& Sigma_cluster_;

  const DualGF& Sigma_dual_;

  DualGF Sigma_lattice_nondiag_;  // \Sigma(\vec{K}, \vec{K'}, \tilde{\vec{k}}, i\omega_n)
  SpLatticeGF& Sigma_lattice_;    // \Sigma(\vec{k} = \vec{K} + \tilde{vec{k}}, i\omega_n)
};

template <typename Scalar, typename Concurrency, int dimension>
void DualToLatticeSelfEnergy<Scalar, Concurrency, dimension>::computeNonDiagonalLatticeSelfEnergy() {
  linalg::Matrix<Complex, linalg::CPU> one_plus_sigma_tilde_G_inv("one_plus_sigma_tilde_G_inv",
                                                                  KClusterDmn::dmn_size());
  linalg::Matrix<Complex, linalg::CPU> sigma_tilde("sigma_tilde", KClusterDmn::dmn_size());
  linalg::Matrix<Complex, linalg::CPU> sigma_bar("sigma_bar", KClusterDmn::dmn_size());

  // Work space for the inverse.
  linalg::Vector<int, linalg::CPU> ipiv;
  linalg::Vector<Complex, linalg::CPU> work;

  // Distribute the work amongst the processes.
  const func::dmn_variadic<KSuperlatticeDmn, SpFreqDmn> k_w_dmn_obj;
  const std::pair<int, int> bounds = concurrency_.get_bounds(k_w_dmn_obj);
  int coor[2];

  Sigma_lattice_nondiag_ = 0.;

  for (int l = bounds.first; l < bounds.second; ++l) {
    k_w_dmn_obj.linind_2_subind(l, coor);
    const auto k_tilde = coor[0];
    const auto w = coor[1];

    for (int K2 = 0; K2 < KClusterDmn::dmn_size(); ++K2) {
      for (int K1 = 0; K1 < KClusterDmn::dmn_size(); ++K1) {
        sigma_tilde(K1, K2) = Sigma_dual_(K1, K2, k_tilde, w);

        one_plus_sigma_tilde_G_inv(K1, K2) =
            Sigma_dual_(K1, K2, k_tilde, w) * G_cluster_(0, 0, 0, 0, K2, w);
        if (K1 == K2)
          one_plus_sigma_tilde_G_inv(K1, K2) += 1.;
      }
    }

    linalg::matrixop::inverse(one_plus_sigma_tilde_G_inv, ipiv, work);

    linalg::matrixop::gemm(one_plus_sigma_tilde_G_inv, sigma_tilde, sigma_bar);

    for (int K2 = 0; K2 < KClusterDmn::dmn_size(); ++K2) {
      for (int K1 = 0; K1 < KClusterDmn::dmn_size(); ++K1) {
        Sigma_lattice_nondiag_(K1, K2, k_tilde, w) = sigma_bar(K1, K2);

        if (K1 == K2)
          Sigma_lattice_nondiag_(K1, K2, k_tilde, w) += Sigma_cluster_(0, 0, 0, 0, K1, w);
      }
    }
  }

  concurrency_.sum(Sigma_lattice_nondiag_);
}

template <typename Scalar, typename Concurrency, int dimension>
void DualToLatticeSelfEnergy<Scalar, Concurrency, dimension>::makeTranslationalInvariant(
    const func::function<Complex, func::dmn_variadic<RClusterDmn, RClusterDmn>>& f_in,
    func::function<Complex, RClusterDmn>& f_out) {
  f_out = 0.;

  for (int R2 = 0; R2 < RClusterDmn::dmn_size(); ++R2) {
    for (int R1 = 0; R1 < RClusterDmn::dmn_size(); ++R1) {
      const int R1_min_R2 = RClusterDmn::parameter_type::subtract(R2, R1);
      f_out(R1_min_R2) += f_in(R1, R2);
    }
  }

  f_out /= RClusterDmn::dmn_size();
}

template <typename Scalar, typename Concurrency, int dimension>
void DualToLatticeSelfEnergy<Scalar, Concurrency, dimension>::updateLatticeSelfEnergy(
    const func::function<Complex, KClusterDmn>& Sigma_lattice_K, const int k_tilde, const int w,
    SpLatticeGF& Sigma_lattice) {
  const auto& k_tilde_vec = KSuperlatticeDmn::get_elements()[k_tilde];

  for (int K = 0; K < KClusterDmn::dmn_size(); ++K) {
    const auto& K_vec = KClusterDmn::get_elements()[K];
    const auto K_plus_k_tilde = math::util::add(K_vec, k_tilde_vec);

    const auto k = phys::domains::cluster_operations::index(
        phys::domains::cluster_operations::find_closest_cluster_vector(
            K_plus_k_tilde, KLatticeDmn::get_elements(),
            KLatticeDmn::parameter_type::get_super_basis_vectors()),
        KLatticeDmn::get_elements(), KLatticeDmn::parameter_type::SHAPE);

    Sigma_lattice(0, 0, 0, 0, k, w) = Sigma_lattice_K(K);
    Sigma_lattice(0, 1, 0, 1, k, w) = Sigma_lattice_K(K);
  }
}

}  // namespace df
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DUAL_FERMION_DUAL_TO_LATTICE_SELF_ENERGY_HPP
