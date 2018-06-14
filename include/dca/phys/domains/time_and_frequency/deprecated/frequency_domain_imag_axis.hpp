// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Frequency domain on the imaginary axis.

#ifndef DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_FREQUENCY_DOMAIN_IMAG_AXIS_HPP
#define DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_FREQUENCY_DOMAIN_IMAG_AXIS_HPP

#include <cmath>
#include <string>
#include <vector>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

class frequency_domain_imag_axis {
public:
  typedef double element_type;

  static int get_size() {
    return get_elements().size();
  }

  static std::string get_name() {
    static std::string name = "frequency-domain-imag-axis";
    return name;
  }

  static std::vector<double>& get_elements() {
    static std::vector<double> vec_elements(0, 0.);
    return vec_elements;
  }

  template <typename Writer>
  static void write(Writer& writer);

  template <typename stream_type>
  static void to_JSON(stream_type& ss);

  template <typename parameters_t>
  static void initialize(parameters_t& parameters);
};

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

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_FREQUENCY_DOMAIN_IMAG_AXIS_HPP
