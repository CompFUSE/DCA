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

#ifndef PHYS_LIBRARY_DOMAINS_QUANTUM_DOMAIN_DCA_ITERATION_DOMAIN_H
#define PHYS_LIBRARY_DOMAINS_QUANTUM_DOMAIN_DCA_ITERATION_DOMAIN_H

#include <string>
#include <vector>

class DCA_iteration_domain {
public:
  typedef int element_type;

public:
  static int& get_size();
  static std::string get_name();

  static std::vector<int>& get_elements();

  template <typename parameters_type>
  static void initialize(parameters_type& parameters);

  template <typename Reader>
  static void read(Reader& reader);

  template <typename Writer>
  static void write(Writer& writer);

  template <class stream_type>
  static void to_JSON(stream_type& ss);

private:
  static std::vector<int>& initialize_elements();
};

int& DCA_iteration_domain::get_size() {
  static int size = 0;
  return size;
}

std::string DCA_iteration_domain::get_name() {
  static std::string name = "DCA-iteration-domain";
  return name;
}

std::vector<int>& DCA_iteration_domain::get_elements() {
  static std::vector<int>& v = initialize_elements();
  return v;
}

template <typename parameters_type>
void DCA_iteration_domain::initialize(parameters_type& parameters) {
  get_size() = parameters.get_DCA_iterations();
}

std::vector<int>& DCA_iteration_domain::initialize_elements() {
  static std::vector<int> v(get_size());

  for (int i = 0; i < get_size(); i++)
    v[i] = i;

  return v;
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

#endif  // PHYS_LIBRARY_DOMAINS_QUANTUM_DOMAIN_DCA_ITERATION_DOMAIN_H
