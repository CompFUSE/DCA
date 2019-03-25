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

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
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
  using BandDmn = func::dmn_0<phys::domains::electron_band_domain>;
  using SpinDmn = func::dmn_0<phys::domains::electron_spin_domain>;

  using SpFreqDmn = func::dmn_0<phys::domains::frequency_domain>;

  using RClusterDmn = typename phys::ClusterDomainAliases<dimension>::RClusterDmn;
  using KClusterDmn = typename phys::ClusterDomainAliases<dimension>::KClusterDmn;
  using KLatticeDmn = typename phys::ClusterDomainAliases<dimension>::KSpHostDmn;
  using KSuperlatticeDmn = typename phys::ClusterDomainAliases<dimension>::KSpSuperlatticeDmn;

  using DualGFDmn = func::dmn_variadic<KClusterDmn, KClusterDmn, KSuperlatticeDmn, SpFreqDmn>;
  using DualGF = func::function<std::complex<Scalar>, DualGFDmn>;

  using RealSpaceDualGFDmn =
      func::dmn_variadic<RClusterDmn, RClusterDmn, KSuperlatticeDmn, SpFreqDmn>;
  using RealSpaceDualGF = func::function<std::complex<Scalar>, RealSpaceDualGFDmn>;

  using SpClusterGFDmn =
      func::dmn_variadic<BandDmn, SpinDmn, BandDmn, SpinDmn, KClusterDmn, SpFreqDmn>;
  using SpClusterGF = func::function<std::complex<Scalar>, SpClusterGFDmn>;

  using SpLatticeGFDmn =
      func::dmn_variadic<BandDmn, SpinDmn, BandDmn, SpinDmn, KLatticeDmn, SpFreqDmn>;
  using SpLatticeGF = func::function<std::complex<Scalar>, SpLatticeGFDmn>;

  DualToLatticeSelfEnergy(const Concurrency& concurrency, const SpClusterGF& G_cluster,
                          const SpClusterGF& Sigma_cluster, const DualGF& Sigma_dual,
                          SpLatticeGF& Sigma_lattice)
      : concurrency_(concurrency),
        G_cluster_(G_cluster),
        Sigma_cluster_(Sigma_cluster),
        Sigma_dual_(Sigma_dual),
        Sigma_lattice_(Sigma_lattice) {}

  // Computes the non-diagonal (in K) lattice self-energy from the dual self-energy.
  void computeNonDiagonalLatticeSelfEnergy() {}

  // Computes the cluster real space Fourier transform of the non-diagonal (in K) lattice self-energy.
  void fourierTransfromToRealSpace() {}

  // Computes the diagonal lattice self-energy from the non-diagonal (in I) real space lattice self-energy.
  void fourierTransfromToMomentumSpace() {}

private:
  const Concurrency& concurrency_;

  const SpClusterGF& G_cluster_;
  const SpClusterGF& Sigma_cluster_;

  const DualGF& Sigma_dual_;

  DualGF Sigma_lattice_nondiag_K_;           // \Sigma(\vec{K}, \vec{K'}, \tilde{\vec{k}})
  RealSpaceDualGF Sigma_lattice_nondiag_R_;  // \Sigma(\vec{I}, \vec{J}, \tilde{\vec{k}})
  SpLatticeGF& Sigma_lattice_;               // \Sigma(\vec{k} = \vec{K} + \tilde{vec{k}})
};

}  // namespace df
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DUAL_FERMION_DUAL_TO_LATTICE_SELF_ENERGY_HPP
