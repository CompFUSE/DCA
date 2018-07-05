// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides the electron band domain.

#ifndef DCA_PHYS_DOMAINS_QUANTUM_ELECTRON_BAND_DOMAIN_HPP
#define DCA_PHYS_DOMAINS_QUANTUM_ELECTRON_BAND_DOMAIN_HPP

#include <cassert>
#include <string>
#include <vector>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

class electron_band_domain {
public:
  struct band_element {
    band_element() : number(-1), flavor(-1), a_vec(0, 0) {}

    int number;
    int flavor;
    std::vector<double> a_vec;
  };

  typedef band_element element_type;

  static int& get_size() {
    static int size = 0;
    return size;
  }

  static std::string get_name() {
    static std::string name = "electron-band-domain";
    return name;
  }

  static std::vector<element_type>& get_elements() {
    static std::vector<element_type> elements(get_size());
    return elements;
  }

  template <typename Writer>
  static void write(Writer& writer);

  template <typename parameters_type>
  static void initialize(parameters_type& parameters, int Nb_bands, std::vector<int> flavors,
                         std::vector<std::vector<double>> a_vecs);
};

template <typename Writer>
void electron_band_domain::write(Writer& writer) {
  writer.open_group(get_name());
  writer.execute(get_elements());
  writer.close_group();
}

template <typename parameters_type>
void electron_band_domain::initialize(parameters_type& /*parameters*/, int NB_BANDS,
                                      std::vector<int> /*flavors*/,
                                      std::vector<std::vector<double>> a_vecs) {
  get_size() = NB_BANDS;

  // assert(NB_BANDS == int(flavors.size()));
  assert(NB_BANDS == int(a_vecs.size()));

  for (size_t i = 0; i < get_elements().size(); ++i) {
    get_elements()[i].number = i;
    get_elements()[i].flavor = i;  // flavors[i];
    get_elements()[i].a_vec = a_vecs[i];
  }
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_QUANTUM_ELECTRON_BAND_DOMAIN_HPP
