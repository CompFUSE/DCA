// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides the Brillouin zone cut domain class that is templated on the number of
// interpolation points.

#ifndef DCA_PHYS_DOMAINS_QUANTUM_BRILLOUIN_ZONE_CUT_DOMAIN_HPP
#define DCA_PHYS_DOMAINS_QUANTUM_BRILLOUIN_ZONE_CUT_DOMAIN_HPP

#include <vector>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <int N = 101>
class brillouin_zone_cut_domain {
public:
  static const int INTERPOLATION_POINTS = N;

  typedef brillouin_zone_cut_domain this_type;
  typedef std::vector<double> element_type;

  static int& get_size() {
    static int size = 0;
    return size;
  }

  static std::vector<element_type>& get_elements() {
    static std::vector<std::vector<double>> v(0);
    return v;
  }

  template <class stream_type>
  static void to_JSON(stream_type& ss);
};

template <int N>
template <class stream_type>
void brillouin_zone_cut_domain<N>::to_JSON(stream_type& ss) {
  ss << "\"brillouin_zone_cut_domain\" : [\n";

  for (int i = 0; i < get_size(); i++) {
    ss << "[ ";
    for (std::size_t z = 0; z < get_elements()[i].size(); z++) {
      if (z == get_elements()[i].size() - 1)
        ss << get_elements()[i][z] << "]";
      else
        ss << get_elements()[i][z] << ", ";
    }

    if (i == int(get_elements().size()) - 1)
      ss << "\n";
    else
      ss << ",\n";
  }
  ss << "]\n";
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_QUANTUM_BRILLOUIN_ZONE_CUT_DOMAIN_HPP
