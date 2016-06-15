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

#ifndef PHYS_LIBRARY_DFT_CONNECTION_VASP_VASP_DOMAINS_VASP_ORBITAL_DOMAIN_HPP
#define PHYS_LIBRARY_DFT_CONNECTION_VASP_VASP_DOMAINS_VASP_ORBITAL_DOMAIN_HPP

#include <string>
#include <vector>

#include "phys_library/domains/cluster/cluster_domain.h"

namespace DFT {
namespace VASP {
// DFT::VASP::

class vasp_orbital_domain {
public:
  typedef int element_type;

public:
  static int& get_size();
  static std::string get_name();

  static std::vector<element_type>& get_elements();

  template <typename Reader>
  static void read(Reader& reader);

  template <typename Writer>
  static void write(Writer& writer);

  template <typename parameters_type>
  static void initialize(parameters_type& parameters);

  template <class stream_type>
  static void to_JSON(stream_type& ss);

private:
  static std::vector<element_type>& initialize_elements();
};

int& vasp_orbital_domain::get_size() {
  static int size = 0;
  return size;
}

std::string vasp_orbital_domain::get_name() {
  static std::string name = "vasp-orbital-domain";
  return name;
}

std::vector<vasp_orbital_domain::element_type>& vasp_orbital_domain::get_elements() {
  static std::vector<element_type> elements(get_size());
  return elements;
}

template <typename Writer>
void vasp_orbital_domain::write(Writer& writer) {
  writer.open_group(get_name());
  writer.execute(get_elements());
  writer.close_group();
}

template <typename parameters_type>
void vasp_orbital_domain::initialize(parameters_type& parameters) {
  get_size() = parameters.get_nb_vasp_orbitals();

  for (size_t i = 0; i < get_elements().size(); ++i) {
    get_elements()[i] = i;
  }
}

}  // VASP
}  // DFT

#endif  // PHYS_LIBRARY_DFT_CONNECTION_VASP_VASP_DOMAINS_VASP_ORBITAL_DOMAIN_HPP
