// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// Description

#ifndef PHYS_LIBRARY_DOMAINS_QUANTUM_DOMAIN_BRILLOUIN_ZONE_CUT_H
#define PHYS_LIBRARY_DOMAINS_QUANTUM_DOMAIN_BRILLOUIN_ZONE_CUT_H

#include <string>
#include <vector>

#include "phys_library/domains/Quantum_domain/Brillouin_zone_cut/brillouin_zone_cut_type.hpp"
#include "phys_library/domains/Quantum_domain/Brillouin_zone_cut/BZC_FERMI_SURFACE_SQUARE_2D_LATTICE.h"
#include "phys_library/domains/Quantum_domain/Brillouin_zone_cut/BZC_SQUARE_2D_LATTICE.h"

template <BrillouinZoneCutType brillouin_zone_cut>
class brillouin_zone_path_domain {
public:
  template <class stream_type>
  static void to_JSON(stream_type& ss);

  template <class stream_type>
  static void to_JSON(stream_type& ss, std::string name, std::vector<std::vector<double>> elements);
};

template <BrillouinZoneCutType brillouin_zone_cut>
template <class stream_type>
void brillouin_zone_path_domain<brillouin_zone_cut>::to_JSON(stream_type& ss) {
  to_JSON(ss, brillouin_zone_path_domain<SQUARE_2D_LATTICE>::get_name(),
          brillouin_zone_path_domain<SQUARE_2D_LATTICE>::get_elements());

  ss << ",\n";
  to_JSON(ss, brillouin_zone_path_domain<FERMI_SURFACE_SQUARE_2D_LATTICE>::get_name(),
          brillouin_zone_path_domain<FERMI_SURFACE_SQUARE_2D_LATTICE>::get_elements());
}

template <BrillouinZoneCutType brillouin_zone_cut>
template <class stream_type>
void brillouin_zone_path_domain<brillouin_zone_cut>::to_JSON(
    stream_type& ss, std::string name, std::vector<std::vector<double>> elements) {
  ss << "\"" << name << "\" : [\n";

  for (size_t i = 0; i < elements.size(); i++) {
    ss << "[ ";
    for (size_t z = 0; z < elements[i].size(); z++) {
      if (z == elements[i].size() - 1)
        ss << elements[i][z] << "]";
      else
        ss << elements[i][z] << ", ";
    }

    if (i == elements.size() - 1)
      ss << "\n";
    else
      ss << ",\n";
  }
  ss << "]\n";
}

#endif  // PHYS_LIBRARY_DOMAINS_QUANTUM_DOMAIN_BRILLOUIN_ZONE_CUT_H
