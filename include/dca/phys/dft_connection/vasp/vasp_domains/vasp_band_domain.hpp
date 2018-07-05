// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// VASP band domain.

#ifndef DCA_PHYS_DFT_CONNECTION_VASP_VASP_DOMAINS_VASP_BAND_DOMAIN_HPP
#define DCA_PHYS_DFT_CONNECTION_VASP_VASP_DOMAINS_VASP_BAND_DOMAIN_HPP

#include <string>
#include <vector>

namespace dca {
namespace phys {
namespace dft {
namespace vasp {
// dca::phys::dft::vasp::

class vasp_band_domain {
public:
  typedef int element_type;

  static int& get_size() {
    static int size = 0;
    return size;
  }

  static std::string get_name() {
    static std::string name = "vasp-band-domain";
    return name;
  }

  static std::vector<element_type>& get_elements() {
    static std::vector<element_type> elements(get_size());
    return elements;
  }

  template <typename Writer>
  static void write(Writer& writer);

  template <typename parameters_type>
  static void initialize(parameters_type& parameters);
};

template <typename Writer>
void vasp_band_domain::write(Writer& writer) {
  writer.open_group(get_name());
  writer.execute(get_elements());
  writer.close_group();
}

template <typename parameters_type>
void vasp_band_domain::initialize(parameters_type& parameters) {
  get_size() = parameters.get_nb_vasp_bands();

  for (size_t i = 0; i < get_elements().size(); ++i) {
    get_elements()[i] = i;
  }
}

}  // vasp
}  // dft
}  // phys
}  // dca

#endif  // DCA_PHYS_DFT_CONNECTION_VASP_VASP_DOMAINS_VASP_BAND_DOMAIN_HPP
