#ifndef DCA_PHYS_DOMAINS_DOMAIN_ALIASES_HPP
#define DCA_PHYS_DOMAINS_DOMAIN_ALIASES_HPP

#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"

namespace dca::phys {
template <class Parameters>
struct DcaDomainAliases {
  using Real = typename Parameters::Real;
  using Scalar = typename Parameters::Scalar;
  using TpAccumulatorPrec = typename Parameters::TPAccumPrec;
  using TpComplex = std::complex<TpAccumulatorPrec>;
  using TDmn = func::dmn_0<domains::time_domain>;
  using WDmn = func::dmn_0<domains::frequency_domain>;
  using WVertexDmn = func::dmn_0<domains::vertex_frequency_domain<domains::COMPACT>>;
  using WExchangeDmn = func::dmn_0<domains::FrequencyExchangeDomain>;

  using BDmn = func::dmn_0<domains::electron_band_domain>;
  using SDmn = func::dmn_0<domains::electron_spin_domain>;
  using NuDmn = func::dmn_variadic<BDmn, SDmn>;  // orbital-spin index
  using NuNuDmn = func::dmn_variadic<NuDmn, NuDmn>;
  static constexpr int Dimension = Parameters::lattice_type::DIMENSION;
  using CDA = ClusterDomainAliases<Dimension>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using KClusterType = typename CDA::KClusterType;
  using KClusterDmn = typename CDA::KClusterDmn;
  using RHostDmn = typename CDA::RSpHostDmn;
  using KHostDmn = typename CDA::KSpHostDmn;
  using KExchangeDmn = func::dmn_0<domains::MomentumExchangeDomain>;
};

}  // namespace dca::phys

#endif
