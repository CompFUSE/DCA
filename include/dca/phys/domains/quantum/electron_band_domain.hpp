// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides the electron band domain.

#ifndef DCA_PHYS_DOMAINS_QUANTUM_ELECTRON_BAND_DOMAIN_HPP
#define DCA_PHYS_DOMAINS_QUANTUM_ELECTRON_BAND_DOMAIN_HPP

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

  static int get_size() {
    return elements_.size();
  }

  static std::string get_name() {
    return "electron-band-domain";
  }

  static const auto& get_elements() {
    return elements_;
  }

  template <typename Writer>
  static void write(Writer& writer);

  template <typename Parameters>
  static void initialize(const Parameters& parameters);

  // For testing purposes only.
  static void setAVectors(const std::vector<std::vector<double>>& vecs){
      if(vecs.size() != get_size()){
          throw(std::logic_error(__PRETTY_FUNCTION__));
      }
      for(int b = 0; b < get_size(); ++b){
          elements_[b].a_vec = vecs[b];
      }
  }

private:
  static inline std::vector<element_type> elements_;
};

template <typename Writer>
void electron_band_domain::write(Writer& writer) {
  writer.open_group(get_name());
  writer.execute(get_elements());
  writer.close_group();
}

template <typename Parameters>
void electron_band_domain::initialize(const Parameters& /*parameters*/) {
  elements_.resize(Parameters::bands);

  using Lattice = typename Parameters::lattice_type;
  auto flavours = Lattice::flavors();
  auto a_vecs = Lattice::aVectors();

  for (size_t i = 0; i < a_vecs.size(); ++i) {
    elements_.at(i).number = i;
    elements_.at(i).flavor = flavours.at(i);
    elements_.at(i).a_vec = a_vecs.at(i);
  }
}

}  // namespace domains
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DOMAINS_QUANTUM_ELECTRON_BAND_DOMAIN_HPP
