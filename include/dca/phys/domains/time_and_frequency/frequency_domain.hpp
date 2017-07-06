// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class parametrizes the fermionic Matsubara frequency domain.
//
// TODO: Use private data members and proper getter methods instead of singletons (see
//       time_domain.hpp).

#ifndef DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_FREQUENCY_DOMAIN_HPP
#define DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_FREQUENCY_DOMAIN_HPP

#include <string>
#include <vector>

#include "dca/math/function_transform/domain_specifications.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

class frequency_domain {
public:
  static constexpr int DIMENSION = 1;

  using scalar_type = double;
  using element_type = double;

  // Needed in function transform.
  using dmn_specifications_type = math::transform::harmonic_dmn_1D_type;

  static bool is_initialized() {
    return initialized_;
  }

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

  // Initializes the elements of the domain with the first num_freqs positive and the first
  // num_freqs negative fermionic Matsubara frequencies using the following order,
  // [-(2*num_freqs-1)*\pi/beta, ..., -\pi/beta, \pi/beta, ..., (2*num_freqs-1)*\pi/beta] .
  static void initialize(double beta, int num_freqs);

  // Calls the previous initialize method with arguments taken from the parameters object.
  template <typename ParametersType>
  static void initialize(const ParametersType& parameters) {
    initialize(parameters.get_beta(), parameters.get_sp_fermionic_frequencies());
  }

private:
  static bool initialized_;
};

template <typename Writer>
void frequency_domain::write(Writer& writer) {
  writer.open_group(get_name());
  writer.execute("elements", get_elements());
  writer.close_group();
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_FREQUENCY_DOMAIN_HPP
