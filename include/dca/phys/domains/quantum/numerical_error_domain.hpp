// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides the numerical error domain.

#ifndef DCA_PHYS_DOMAINS_QUANTUM_NUMERICAL_ERROR_DOMAIN_HPP
#define DCA_PHYS_DOMAINS_QUANTUM_NUMERICAL_ERROR_DOMAIN_HPP

#include <string>
#include <vector>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

class numerical_error_domain {
public:
  typedef double element_type;
  typedef numerical_error_domain this_type;

  static std::string get_name() {
    return "numerical-error-domain";
  }

  static int get_size() {
    return int(get_elements().size());
  }

  static std::vector<element_type>& get_elements() {
    static std::vector<element_type>& v = initialize_elements();
    return v;
  }

  template <typename Writer>
  static void write(Writer& writer);

private:
  static std::vector<element_type>& initialize_elements();
};

template <typename Writer>
void numerical_error_domain::write(Writer& writer) {
  writer.open_group(get_name());
  writer.execute("elements", get_elements());
  writer.close_group();
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_QUANTUM_NUMERICAL_ERROR_DOMAIN_HPP
