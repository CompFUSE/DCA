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
    return elements_;
  }

  template <typename Writer>
  static void write(Writer& writer);

  static int to_coordinate(element_type spin);

  template <typename Parameters>
  static void initialize(const Parameters& parameters);

private:
  static std::vector<e_spin_states_type> initialize_elements();
  static void initialize_elements(std::vector<e_spin_states_type>& elements);
  static inline std::vector<element_type> elements_;
};

template <typename Writer>
void electron_spin_domain::write(Writer& writer) {
  writer.open_group(get_name());
  writer.execute("elements", get_elements());
  writer.close_group();
}

template <typename Parameters>
void electron_spin_domain::initialize(const Parameters& /*parameters*/) {
  using Lattice = typename Parameters::lattice_type;
  int spins = 2;
  if constexpr(Parameters::template HasCustomSpin<Parameters::Model::lattice_type::SPINS>::value) {
    spins = Parameters::Model::lattice_type::SPINS;
  }
  elements_.resize(spins);
  initialize_elements(elements_);
}

  
}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_QUANTUM_ELECTRON_SPIN_DOMAIN_HPP
