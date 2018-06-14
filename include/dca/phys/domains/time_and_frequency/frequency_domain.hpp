// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class parametrizes the fermionic Matsubara frequency domain.

#ifndef DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_FREQUENCY_DOMAIN_HPP
#define DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_FREQUENCY_DOMAIN_HPP

#include <cassert>
#include <cstdlib>  // std::size_t
#include <string>
#include <vector>

#include "dca/math/function_transform/domain_specifications.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

class frequency_domain {
public:
  static constexpr int dimension = 1;

  using ScalarType = double;
  using element_type = ScalarType;

  // Needed in function transform.
  using dmn_specifications_type = math::transform::harmonic_dmn_1D_type;

  static bool is_initialized() {
    return initialized_;
  }

  static const std::string& get_name() {
    return name_;
  }

  static std::size_t get_size() {
    assert(initialized_);
    return elements_.size();
  }

  // TODO: Add const qualifier when rest of the code is fixed.
  static /*const*/ std::vector<element_type>& get_elements() {
    assert(initialized_);
    return elements_;
  }

  // Returns the Matsubara frequency indices of the elements.
  static const std::vector<int>& get_indices() {
    assert(initialized_);
    return indices_;
  }

  template <typename Writer>
  static void write(Writer& writer);

  // Initializes the elements of the domain with the first num_freqs positive and the first
  // num_freqs negative fermionic Matsubara frequencies. The elements are sorted w.r.t the Matsubara
  // frequency index in increasing order,
  // [-(2*num_freqs-1)*\pi/beta, -(2*num_freqs-3)*\pi/beta, ..., -\pi/beta, \pi/beta, ...,
  // (2*num_freqs-1)*\pi/beta] .
  static void initialize(ScalarType beta, int num_freqs);

  // Calls the previous initialize method with arguments taken from the parameters object.
  template <typename ParametersType>
  static void initialize(const ParametersType& parameters) {
    initialize(parameters.get_beta(), parameters.get_sp_fermionic_frequencies());
  }

private:
  static bool initialized_;
  const static std::string name_;
  static std::vector<element_type> elements_;
  static std::vector<int> indices_;
};

template <typename Writer>
void frequency_domain::write(Writer& writer) {
  writer.open_group(name_);
  writer.execute("elements", elements_);
  writer.execute("indices", indices_);
  writer.close_group();
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_FREQUENCY_DOMAIN_HPP
