// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class converts pairs of electron band and electron spin indices into single spin orbital
// indices and vice versa.

#ifndef DCA_PHYS_DOMAINS_CONVERT_HPP
#define DCA_PHYS_DOMAINS_CONVERT_HPP

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

// Empty class template.
template <typename target, typename whatever>
class convert {};

// Specialization for conversion of an electron band index and an electron spin state into a single
// spin orbital index.
template <>
class convert<int, func::dmn_variadic<func::dmn_0<domains::electron_band_domain>,
                                      func::dmn_0<domains::electron_spin_domain>>> {
public:
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin domain

  static int spin_orbital(int band, e_spin_states_type e_spin);

private:
  static func::function<int, nu>& intitialize_spin_orbital();
};

// Specialization for conversion of a single spin orbital index into an electron band and electron
// spin index.
template <>
class convert<func::dmn_variadic<func::dmn_0<domains::electron_band_domain>,
                                 func::dmn_0<domains::electron_spin_domain>>,
              int> {
public:
  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;

  static std::pair<int, int> from_spin_orbital(int spo);
};

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_CONVERT_HPP
