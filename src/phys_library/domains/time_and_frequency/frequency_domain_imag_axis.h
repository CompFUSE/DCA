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

#ifndef PHYS_LIBRARY_DOMAINS_TIME_AND_FREQUENCY_FREQUENCY_DOMAIN_IMAG_AXIS_H
#define PHYS_LIBRARY_DOMAINS_TIME_AND_FREQUENCY_FREQUENCY_DOMAIN_IMAG_AXIS_H

#include <cmath>
#include <string>
#include <vector>

class frequency_domain_imag_axis {
public:
  typedef double element_type;

public:
  static int get_size();
  static std::string get_name();

  static std::vector<double>& get_elements();

  template <typename stream_type>
  static void to_JSON(stream_type& ss);

  template <typename Writer>
  static void write(Writer& writer);

  template <typename parameters_t>
  static void initialize(parameters_t& parameters);
};

int frequency_domain_imag_axis::get_size() {
  return get_elements().size();
}

std::string frequency_domain_imag_axis::get_name() {
  static std::string name = "frequency-domain-imag-axis";
  return name;
}

std::vector<double>& frequency_domain_imag_axis::get_elements() {
  static std::vector<double> vec_elements(0, 0.);
  return vec_elements;
}

template <typename Writer>
void frequency_domain_imag_axis::write(Writer& writer) {
  writer.open_group(get_name());

  writer.execute("elements", get_elements());

  writer.close_group();
}

template <typename stream_type>
void frequency_domain_imag_axis::to_JSON(stream_type& ss) {
  ss << "\"Pade_imag_frequency_domain\" : [\n";

  for (int i = 0; i < get_size(); i++)
    if (i == get_size() - 1)
      ss << get_elements()[i] << "\n";
    else
      ss << get_elements()[i] << ",\n";
  ss << "]\n";
}

template <typename parameters_t>
void frequency_domain_imag_axis::initialize(parameters_t& parameters) {
  get_elements().resize(parameters.get_N_wn());

  for (int l = 0; l < get_size(); l++)
    get_elements()[l] = M_PI / parameters.get_beta() * (1 + 2 * l);
}

#endif  // PHYS_LIBRARY_DOMAINS_TIME_AND_FREQUENCY_FREQUENCY_DOMAIN_IMAG_AXIS_H
