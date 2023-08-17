// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides the electron spin domain.

#ifndef DCA_PHYS_DOMAINS_QUANTUM_ELECTRON_SPIN_DOMAIN_HPP
#define DCA_PHYS_DOMAINS_QUANTUM_ELECTRON_SPIN_DOMAIN_HPP

#include <string>
#include <vector>

#include "dca/phys/domains/quantum/e_spin_states.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

class electron_spin_domain {
public:
  typedef e_spin_states_type element_type;

  static int get_size() {
    return 2;
  }

  static std::string get_name() {
    static std::string name = "electron-spin-domain";
    return name;
  }

  static std::vector<e_spin_states_type>& get_elements() {
    static std::vector<e_spin_states_type> v = initialize_elements();
    return v;
  }

  template <typename Writer>
  static void write(Writer& writer);

  static int to_coordinate(element_type spin);

private:
  static std::vector<e_spin_states_type> initialize_elements();
};

template <typename Writer>
void electron_spin_domain::write(Writer& writer) {
  writer.open_group(get_name());
  writer.execute("elements", get_elements());
  writer.close_group();
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_QUANTUM_ELECTRON_SPIN_DOMAIN_HPP
