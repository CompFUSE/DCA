// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides the DCA iteration domain.

#ifndef DCA_PHYS_DOMAINS_QUANTUM_DCA_ITERATION_DOMAIN_HPP
#define DCA_PHYS_DOMAINS_QUANTUM_DCA_ITERATION_DOMAIN_HPP

#include <string>
#include <vector>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

class DCA_iteration_domain {
public:
  typedef int element_type;

  static int& get_size() {
    static int size = 0;
    return size;
  }

  static std::string get_name() {
    static std::string name = "DCA-iteration-domain";
    return name;
  }

  static std::vector<int>& get_elements() {
    static std::vector<int>& v = initialize_elements();
    return v;
  }

  template <typename parameters_type>
  static void initialize(parameters_type& parameters);

  template <typename Writer>
  static void write(Writer& writer);

  template <class stream_type>
  static void to_JSON(stream_type& ss);

private:
  static std::vector<int>& initialize_elements();
};

template <typename parameters_type>
void DCA_iteration_domain::initialize(parameters_type& parameters) {
  get_size() = parameters.get_dca_iterations();
}

template <typename Writer>
void DCA_iteration_domain::write(Writer& writer) {
  writer.open_group(get_name());
  writer.execute("elements", get_elements());
  writer.close_group();
}

template <class stream_type>
void DCA_iteration_domain::to_JSON(stream_type& ss) {
  ss << "\"DCA_iteration_domain\" : [\n";

  for (int i = 0; i < get_size(); i++)
    if (i == get_size() - 1)
      ss << get_elements()[i] << "\n";
    else
      ss << get_elements()[i] << ",\n";

  ss << "]\n";
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_QUANTUM_DCA_ITERATION_DOMAIN_HPP
