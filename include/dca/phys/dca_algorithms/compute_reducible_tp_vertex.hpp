// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides a method to compute the reducible two-particle vertex in the channels
// - particle-hole longitudinal up-up,
// - particle-hole longitudinal up-down, and
// - particle-hole transverse.

#ifndef DCA_PHYS_DCA_ALGORITHMS_COMPUTE_REDUCIBLE_TP_VERTEX_HPP
#define DCA_PHYS_DCA_ALGORITHMS_COMPUTE_REDUCIBLE_TP_VERTEX_HPP

#include <cassert>
#include <complex>
#include <stdexcept>
#include <utility>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/four_point_type.hpp"

namespace dca {
namespace phys {
// dca::phys::

template <typename Scalar, typename OrbitalDmn, typename SpinDmn, typename SpFreqDmn,
          typename TpFreqDmn, typename FreqExchangeDmn, typename KDmn, typename ConcurrencyType>
void computeReducibleTpVertex(
    const Scalar beta, const ConcurrencyType& concurrency,
    const func::function<std::complex<Scalar>,
                         func::dmn_variadic<OrbitalDmn, SpinDmn, OrbitalDmn, SpinDmn, KDmn, SpFreqDmn>>& G,
    const func::function<std::complex<Scalar>,
                         func::dmn_variadic<OrbitalDmn, OrbitalDmn, OrbitalDmn, OrbitalDmn, KDmn,
                                            KDmn, KDmn, TpFreqDmn, TpFreqDmn, FreqExchangeDmn>>& G4,
    const FourPointType channel,
    func::function<std::complex<Scalar>,
                   func::dmn_variadic<OrbitalDmn, OrbitalDmn, OrbitalDmn, OrbitalDmn, KDmn, KDmn,
                                      KDmn, TpFreqDmn, TpFreqDmn, FreqExchangeDmn>>& Gamma) {
  // TODO: Multi-orbital support.
  assert(OrbitalDmn::dmn_size() == 1);

  if (!(channel == PARTICLE_HOLE_LONGITUDINAL_UP_UP ||
        channel == PARTICLE_HOLE_LONGITUDINAL_UP_DOWN || channel == PARTICLE_HOLE_TRANSVERSE))
    throw std::invalid_argument(
        "The computation of the reducible two-particle vertex is only implemented for the channels "
        "particle-hole longitudinal up-up, longitudinal up-down, and transverse.");

  auto k_plus_q = [](const int k, const int q) { return KDmn::parameter_type::add(k, q); };
  auto w_plus_w_ex = [](const int w, const int w_ex) { return w + w_ex; };

  const auto Nc = KDmn::dmn_size();
  const auto k0 = KDmn::parameter_type::origin_index();
  const auto w_ex_0 = 0;

  // Need to set everything to zero in the beginning for the final reduction to work.
  Gamma = 0.;

  // Distribute the work amongst the processes.
  const func::dmn_variadic<KDmn, FreqExchangeDmn> k_w_ex_dmn_obj;
  const std::pair<int, int> bounds = concurrency.get_bounds(k_w_ex_dmn_obj);
  int coor[2];

  for (int l = bounds.first; l < bounds.second; ++l) {
    k_w_ex_dmn_obj.linind_2_subind(l, coor);
    const auto q = coor[0];
    const auto w_ex = coor[1];

    for (int w2 = 0; w2 < TpFreqDmn::dmn_size(); ++w2) {
      const auto w2_plus_w_ex = w_plus_w_ex(w2, w_ex);
      for (int w1 = 0; w1 < TpFreqDmn::dmn_size(); ++w1) {
        const auto w1_plus_w_ex = w_plus_w_ex(w1, w_ex);

        for (int k2 = 0; k2 < KDmn::dmn_size(); ++k2) {
          const auto k2_plus_q = k_plus_q(k2, q);
          for (int k1 = 0; k1 < KDmn::dmn_size(); ++k1) {
            const auto k1_plus_q = k_plus_q(k1, q);

            Gamma(0, 0, 0, 0, k1, k2, q, w1, w2, w_ex) +=
                beta * Nc * G4(0, 0, 0, 0, k1, k2, q, w1, w2, w_ex) /
                (G(0, 0, 0, 0, k1_plus_q, w1_plus_w_ex) * G(0, 0, 0, 0, k1, w1) *
                 G(0, 0, 0, 0, k2_plus_q, w2_plus_w_ex) * G(0, 0, 0, 0, k2, w2));

            // Vertical legs.
            if (channel == PARTICLE_HOLE_LONGITUDINAL_UP_UP ||
                channel == PARTICLE_HOLE_LONGITUDINAL_UP_DOWN) {
              if (q == k0 && w_ex == w_ex_0) {
                Gamma(0, 0, 0, 0, k1, k2, q, w1, w2, w_ex) -=
                    beta * Nc / (G(0, 0, 0, 0, k1, w1) * G(0, 0, 0, 0, k2, w2));
              }
            }

            // Horizontal legs.
            if (channel == PARTICLE_HOLE_LONGITUDINAL_UP_UP || channel == PARTICLE_HOLE_TRANSVERSE) {
              if (k1 == k2 && w1 == w2) {
                Gamma(0, 0, 0, 0, k1, k2, q, w1, w2, w_ex) +=
                    beta * Nc / (G(0, 0, 0, 0, k1_plus_q, w1_plus_w_ex) * G(0, 0, 0, 0, k1, w1));
              }
            }
          }
        }
      }
    }
  }

  concurrency.sum(Gamma);
}

}  // phys
}  // dca

#endif  //  DCA_PHYS_DCA_ALGORITHMS_COMPUTE_REDUCIBLE_TP_VERTEX_HPP
