// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Frequency domain on the real axis.

#ifndef DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_FREQUENCY_DOMAIN_REAL_AXIS_HPP
#define DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_FREQUENCY_DOMAIN_REAL_AXIS_HPP

#include <string>
#include <vector>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

class frequency_domain_real_axis {
public:
  typedef double element_type;

  static int get_size() {
    return get_elements().size();
  }

  static std::string get_name() {
    static std::string name = "frequency-domain-real-axis";
    return name;
  }

  static std::vector<double>& get_elements() {
    static std::vector<double> vec_elements(0, 0);
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
void frequency_domain_real_axis::write(Writer& writer) {
  writer.open_group(get_name());
  writer.execute("elements", get_elements());
  writer.close_group();
}

template <typename stream_type>
void frequency_domain_real_axis::to_JSON(stream_type& ss) {
  ss << "\"Pade_real_frequency_domain\" : [\n";

  for (int i = 0; i < get_size(); i++)
    if (i == get_size() - 1)
      ss << get_elements()[i] << "\n";
    else
      ss << get_elements()[i] << ",\n";
  ss << "]\n";
}

template <typename parameters_t>
void frequency_domain_real_axis::initialize(parameters_t& parameters) {
  int N = parameters.get_real_frequencies();

  double min = parameters.get_min_real_frequency();
  double max = parameters.get_max_real_frequency();

  double delta = (max - min) / double(N - 1);

  get_elements().resize(N, min);

  for (int l = 0; l < get_size(); l++)
    get_elements()[l] += double(l) * delta;
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_FREQUENCY_DOMAIN_REAL_AXIS_HPP
