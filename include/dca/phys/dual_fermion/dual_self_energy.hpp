// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class implements the dual self-energy up to second order.

#ifndef DCA_PHYS_DUAL_FERMION_DUAL_SELF_ENERGY_HPP
#define DCA_PHYS_DUAL_FERMION_DUAL_SELF_ENERGY_HPP

#include <cassert>
#include <complex>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"

namespace dca {
namespace phys {
namespace df {
// dca::phys::df::

template <typename Scalar, typename Concurrency, typename BandDmn, typename KClusterDmn,
          typename KSuperlatticeDmn, typename TpFreqDmn, typename FreqExchangeDmn>
class DualSelfEnergy {
public:
  using TpGreensFunctionDomain =
      func::dmn_variadic<BandDmn, BandDmn, BandDmn, BandDmn, KClusterDmn, KClusterDmn, KClusterDmn,
                         TpFreqDmn, TpFreqDmn, FreqExchangeDmn>;
  using TpGreensFunction = func::function<std::complex<Scalar>, TpGreensFunctionDomain>;

  using DualGreensFunctionDomain =
      func::dmn_variadic<KClusterDmn, KClusterDmn, KSuperlatticeDmn, TpFreqDmn>;
  using DualGreensFunction = func::function<std::complex<Scalar>, DualGreensFunctionDomain>;

  DualSelfEnergy(const Concurrency& concurrency, const Scalar beta,
                 const DualGreensFunction& G0_tilde, const TpGreensFunction& Gamma_uu,
                 const TpGreensFunction& Gamma_ud)
      : concurrency_(concurrency),
        beta_(beta),
        G0_tilde_(G0_tilde),
        Gamma_uu_(Gamma_uu),
        Gamma_ud_(Gamma_ud) {
    // TODO: Multi-orbital support.
    assert(BandDmn::dmn_size() == 1);
  }

  void compute1stOrder(){};

  const DualGreensFunction& get() {
    return Sigma_tilde_;
  }

private:
  const Concurrency& concurrency_;
  const Scalar beta_;

  const DualGreensFunction& G0_tilde_;
  const TpGreensFunction& Gamma_uu_;
  const TpGreensFunction& Gamma_ud_;

  DualGreensFunction Sigma_tilde_;
};

}  // namespace df
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DUAL_FERMION_DUAL_SELF_ENERGY_HPP
