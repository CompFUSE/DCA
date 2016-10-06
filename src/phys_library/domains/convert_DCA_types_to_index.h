// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements various types and orderings of the matsubara-frequency domain via
// templates.

#ifndef PHYS_LIBRARY_DOMAINS_CONVERT_DCA_TYPES_TO_INDEX_H
#define PHYS_LIBRARY_DOMAINS_CONVERT_DCA_TYPES_TO_INDEX_H

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"

namespace QMC {
template <typename target, class whatever>
class convert {};
template <>
class convert<
    int, func::dmn_variadic<func::dmn_0<electron_band_domain>, func::dmn_0<electron_spin_domain>>> {
public:
  using b = func::dmn_0<electron_band_domain>;
  using s = func::dmn_0<electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin domain

  static int spin_orbital(int band, e_spin_states_type e_spin);

private:
  static func::function<int, nu>& intitialize_spin_orbital();
};

int convert<int, func::dmn_variadic<func::dmn_0<electron_band_domain>,
                                    func::dmn_0<electron_spin_domain>>>::spin_orbital(int band,
                                                                                      e_spin_states_type
                                                                                          e_spin) {
  typedef nu::this_type parameter_typelist;

  static func::function<int, nu>& spo_function = intitialize_spin_orbital();

  if (std::is_same<dca::util::TypeAt<0, parameter_typelist>::type, electron_band_domain>::value) {
    return spo_function(band, electron_spin_domain::to_coordinate(e_spin));
  }
  else {
    return spo_function(electron_spin_domain::to_coordinate(e_spin), band);
  }
}

func::function<int, func::dmn_variadic<func::dmn_0<electron_band_domain>, func::dmn_0<electron_spin_domain>>>& convert<
    int, func::dmn_variadic<func::dmn_0<electron_band_domain>,
                            func::dmn_0<electron_spin_domain>>>::intitialize_spin_orbital() {
  static func::function<int, nu> spo_function;

  for (int i = 0; i < spo_function.size(); i++)
    spo_function(i) = i;

  return spo_function;
}

template <>
class convert<func::dmn_variadic<func::dmn_0<electron_band_domain>, func::dmn_0<electron_spin_domain>>,
              int> {
public:
  using b = func::dmn_0<electron_band_domain>;
  using s = func::dmn_0<electron_spin_domain>;

  static std::pair<int, int> from_spin_orbital(int spo);
};

std::pair<int, int> convert<
    func::dmn_variadic<func::dmn_0<electron_band_domain>, func::dmn_0<electron_spin_domain>>,
    int>::from_spin_orbital(int spo) {
  std::pair<int, int> tmp;
  tmp.first = spo % b::dmn_size();
  tmp.second = (spo - tmp.first) % s::dmn_size();
  return tmp;
}
}

#endif  // PHYS_LIBRARY_DOMAINS_CONVERT_DCA_TYPES_TO_INDEX_H
