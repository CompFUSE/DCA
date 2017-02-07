// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Frequency domain.

#ifndef DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_FREQUENCY_DOMAIN_HPP
#define DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_FREQUENCY_DOMAIN_HPP

#include <cmath>
#include <string>
#include <vector>

#include "dca/math/function_transform/domain_specifications.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

class frequency_domain {
public:
  const static int DIMENSION = 1;

  typedef double scalar_type;
  typedef double element_type;

  typedef math::transform::harmonic_dmn_1D_type dmn_specifications_type;

  static int& get_size() {
    static int size;
    return size;
  }

  static std::string get_name() {
    static std::string name = "frequency-domain";
    return name;
  }

  static scalar_type* get_basis() {
    static scalar_type basis[DIMENSION];
    return basis;
  }

  static scalar_type* get_inverse_basis() {
    static scalar_type inv_basis[DIMENSION];
    return inv_basis;
  }

  static std::vector<double>& get_elements() {
    static std::vector<double> elements;
    return elements;
  }

  static std::vector<int>& get_integer_wave_vectors() {
    static std::vector<int> elements;
    return elements;
  }

  template <typename Writer>
  static void write(Writer& writer);

  static void initialize(double beta, int n_frequencies);

  template <typename parameters_t>
  static void initialize(parameters_t& parameters){
    initialize(parameters.get_beta(), parameters.get_sp_fermionic_frequencies());
  }
};

template <typename Writer>
void frequency_domain::write(Writer& writer) {
  writer.open_group(get_name());
  writer.execute("elements", get_elements());
  writer.close_group();
}

void frequency_domain::initialize(double beta, int n_frequencies) {
  get_basis()[0] = (2. * M_PI) / beta;
  get_inverse_basis()[0] = beta / (2. * M_PI);

  get_size() = 2 * n_frequencies;

  get_elements().resize(get_size());
  get_integer_wave_vectors().resize(get_size());

  for (int l = 0; l < n_frequencies; l++) {
    get_elements()[get_size() / 2 + 0 + l] = M_PI / beta * (1 + 2 * l);
    get_elements()[get_size() / 2 - 1 - l] = -M_PI / beta * (1 + 2 * l);

    get_integer_wave_vectors()[get_size() / 2 + 0 + l] = (1 + 2 * l);
    get_integer_wave_vectors()[get_size() / 2 - 1 - l] = -(1 + 2 * l);
  }
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_FREQUENCY_DOMAIN_HPP
