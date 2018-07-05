// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides the point group symmetry domain.

#ifndef DCA_PHYS_DOMAINS_QUANTUM_POINT_GROUP_SYMMETRY_DOMAIN_HPP
#define DCA_PHYS_DOMAINS_QUANTUM_POINT_GROUP_SYMMETRY_DOMAIN_HPP

#include "dca/phys/domains/quantum/point_group_symmetry_element.hpp"
#include "dca/phys/domains/quantum/symmetry_group_level.hpp"

#include <cassert>
#include <string>
#include <vector>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <symmetry_group_level_type symmetry_group_level, typename base_cluster_type>
class point_group_symmetry_domain {
public:
  typedef point_group_symmetry_element element_type;

  static int DIMENSION;

  static int& get_size() {
    static int size = 0;
    return size;
  }

  static std::string get_name() {
    static std::string name = "point-group-symmetry-domain (" + to_str(symmetry_group_level) +
                              ", " + base_cluster_type::get_name() + ")";
    return name;
  }

  static std::vector<element_type>& get_elements() {
    static std::vector<point_group_symmetry_element> v(get_size(),
                                                       point_group_symmetry_element(DIMENSION));
    return v;
  }

  static bool is_initialized() {
    static bool initialized = false;
    return initialized;
  }

  static std::vector<double> transform(int i, std::vector<double> vec);

  template <typename Writer>
  static void write(Writer& writer);

  template <class stream_type>
  static void to_JSON(stream_type& ss);
};

template <symmetry_group_level_type symmetry_group_level, typename base_cluster_type>
int point_group_symmetry_domain<symmetry_group_level, base_cluster_type>::DIMENSION = -1;

template <symmetry_group_level_type symmetry_group_level, typename base_cluster_type>
std::vector<double> point_group_symmetry_domain<symmetry_group_level, base_cluster_type>::transform(
    int i, std::vector<double> vec) {
  assert(i > -1);
  assert(i < get_size());

  std::vector<double> tmp(DIMENSION, 0);

  get_elements()[i].transform(&vec[0], &tmp[0]);

  return tmp;
}

template <symmetry_group_level_type symmetry_group_level, typename base_cluster_type>
template <typename Writer>
void point_group_symmetry_domain<symmetry_group_level, base_cluster_type>::write(Writer& /*writer*/) {
  //   std::vector<double>               T(DIMENSION);
  //   std::vector<std::vector<double> > O(DIMENSION, std::vector<double>(DIMENSION));

  //   writer.open_group(get_name());

  //   writer.close_group();
}

template <symmetry_group_level_type symmetry_group_level, typename base_cluster_type>
template <class stream_type>
void point_group_symmetry_domain<symmetry_group_level, base_cluster_type>::to_JSON(stream_type& ss) {
  ss << "\"pointgroup symmetry domain ";

  if (symmetry_group_level == UNIT_CELL)
    ss << "unit-cell";

  if (symmetry_group_level == SUPER_CELL)
    ss << "super-cell";

  ss << "\" : {\n";

  for (int l = 0; l < get_size(); ++l) {
    ss << "\"" << l << "\" : \n{\n";

    get_elements()[l].to_JSON(ss);

    if (l == get_size() - 1)
      ss << "\n}\n";
    else
      ss << "\n},\n";
  }
  ss << "}\n";
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_QUANTUM_POINT_GROUP_SYMMETRY_DOMAIN_HPP
