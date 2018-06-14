// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides the particle number domain.

#ifndef DCA_PHYS_DOMAINS_QUANTUM_PARTICLE_NUMBER_DOMAIN_HPP
#define DCA_PHYS_DOMAINS_QUANTUM_PARTICLE_NUMBER_DOMAIN_HPP

#include <vector>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <typename band_dmn_t, typename cluster_t, typename e_spin_t>
class particle_number_domain {
public:
  typedef int element_type;

  static int& get_size() {
    static int size = band_dmn_t::dmn_size() * cluster_t::dmn_size() * e_spin_t::dmn_size() + 1;
    return size;
  }

  static std::vector<int>& get_elements() {
    static std::vector<int>& v = initialize_elements();
    return v;
  }

private:
  static std::vector<int>& initialize_elements();
};

template <typename band_dmn_t, typename cluster_t, typename e_spin_t>
std::vector<int>& particle_number_domain<band_dmn_t, cluster_t, e_spin_t>::initialize_elements() {
  static std::vector<int> v(get_size());

  for (int i = 0; i < get_size(); i++)
    v[i] = i;

  return v;
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_QUANTUM_PARTICLE_NUMBER_DOMAIN_HPP
