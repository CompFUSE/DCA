#ifndef DCA_PHYS_TYPES_SHARED_TYPES_HPP
#define DCA_PHYS_TYPES_SHARED_TYPES_HPP

#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/phys/domains/domain_aliases.hpp"
#include "dca/function/function.hpp"
#include "dca/function/domains/dmn_variadic.hpp"

namespace dca::phys {
template <class Parameters>
struct DcaSharedTypes {
  using DDA = DcaDomainAliases<Parameters>;
  using Real = typename DDA::Real;
  using DisorderConfiguration =
      func::function<Real, func::dmn_variadic<typename DDA::NuDmn, typename DDA::RClusterDmn>>;
};
}  // namespace dca::phys
#endif
