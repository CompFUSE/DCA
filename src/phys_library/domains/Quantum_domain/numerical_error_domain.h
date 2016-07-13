// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef PHYS_LIBRARY_DOMAINS_QUANTUM_DOMAIN_NUMERICAL_ERROR_DOMAIN_H
#define PHYS_LIBRARY_DOMAINS_QUANTUM_DOMAIN_NUMERICAL_ERROR_DOMAIN_H

#include <cmath>
#include <string>
#include <vector>

class numerical_error_domain {
public:
  typedef double element_type;
  typedef numerical_error_domain this_type;

public:
  static std::string get_name();

  static int get_size();
  static std::vector<element_type>& get_elements();

  template <typename Writer>
  static void write(Writer& writer);

private:
  static std::vector<element_type>& initialize_elements();
};

std::string numerical_error_domain::get_name() {
  return "numerical-error-domain";
}

int numerical_error_domain::get_size() {
  return int(get_elements().size());
}

std::vector<double>& numerical_error_domain::get_elements() {
  static std::vector<element_type>& v = initialize_elements();
  return v;
}

std::vector<double>& numerical_error_domain::initialize_elements() {
  static std::vector<element_type> v(0);

  for (int i = -16; i < 0; i++)
    for (double j = 1; j < 10; j += 1.)
      v.push_back(j * std::pow(10., i));

  return v;
}

template <typename Writer>
void numerical_error_domain::write(Writer& writer) {
  writer.open_group(get_name());
  writer.execute("elements", get_elements());
  writer.close_group();
}

#endif  // PHYS_LIBRARY_DOMAINS_QUANTUM_DOMAIN_NUMERICAL_ERROR_DOMAIN_H
