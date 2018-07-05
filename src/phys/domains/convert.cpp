// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements convert.hpp.

#include "dca/phys/domains/convert.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

int convert<int, func::dmn_variadic<func::dmn_0<domains::electron_band_domain>,
                                    func::dmn_0<domains::electron_spin_domain>>>::
    spin_orbital(int band, e_spin_states_type e_spin) {
  typedef nu::this_type parameter_typelist;

  static func::function<int, nu>& spo_function = intitialize_spin_orbital();

  if (std::is_same<dca::util::TypeAt<0, parameter_typelist>::type, domains::electron_band_domain>::value) {
    return spo_function(band, domains::electron_spin_domain::to_coordinate(e_spin));
  }
  else {
    return spo_function(domains::electron_spin_domain::to_coordinate(e_spin), band);
  }
}

func::function<int, func::dmn_variadic<func::dmn_0<domains::electron_band_domain>, func::dmn_0<domains::electron_spin_domain>>>& convert<
    int, func::dmn_variadic<func::dmn_0<domains::electron_band_domain>,
                            func::dmn_0<domains::electron_spin_domain>>>::intitialize_spin_orbital() {
  static func::function<int, nu> spo_function;

  for (int i = 0; i < spo_function.size(); i++)
    spo_function(i) = i;

  return spo_function;
}

std::pair<int, int> convert<func::dmn_variadic<func::dmn_0<domains::electron_band_domain>,
                                               func::dmn_0<domains::electron_spin_domain>>,
                            int>::from_spin_orbital(int spo) {
  std::pair<int, int> tmp;
  tmp.first = spo % b::dmn_size();
  tmp.second = (spo - tmp.first) % s::dmn_size();
  return tmp;
}

}  // domains
}  // phys
}  // dca
