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
// This class parametrizes the imaginary time domain with elements defined in the interval
// [-beta, beta].
// See the description of the first initialize member function for how the elements are initialized.

#ifndef DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_TIME_DOMAIN_HPP
#define DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_TIME_DOMAIN_HPP

#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>

#include "dca/math/function_transform/domain_specifications.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

class time_domain {
public:
  using element_type = double;
  using scalar_type = double;

  // Needed in function transform.
  using dmn_specifications_type = math::transform::interval_dmn_1D_type;

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

  // INTERNAL: Why is the volume = beta?
  static double get_volume() {
    assert(initialized_);
    return beta_;
  }

  template <typename Writer>
  static void write(Writer& writer);

  template <class stream_type>
  static void to_JSON(stream_type& ss);

  // Initializes the elements of the imaginary time domain as follows,
  // [-beta+eps, -beta+step, ..., -step, -eps, eps, step, ..., beta-step, beta-eps],
  // where step = beta/time_slices.
  // Note that the boundaries and the two zeros are shifted by a small number 'eps'.
  static void initialize(const scalar_type beta, const int time_slices,
                         const scalar_type eps = 1.e-16);

  // Calls the previous initialize method with arguments taken from the parameters object.
  template <typename parameters_t>
  static void initialize(const parameters_t& parameters) {
    initialize(parameters.get_beta(), parameters.get_sp_time_intervals());
  }

  // Initializes weights and nodes for integration on the imaginary time domain.
  // The number of nodes per time slice is 2^level.
  // Precondition: The time domain is initialized, i.e. one of the initialize methods has been
  //               called.
  // INTERNAL: Is this method deprecated?
  static void initialize_integration_domain(const int level, std::vector<scalar_type>& weights,
                                            std::vector<element_type>& nodes);

private:
  static bool initialized_;
  const static std::string name_;
  static scalar_type beta_;
  static std::vector<element_type> elements_;
};

template <typename Writer>
void time_domain::write(Writer& writer) {
  writer.open_group(name_);
  writer.execute("elements", elements_);
  writer.close_group();
}

template <class stream_type>
void time_domain::to_JSON(stream_type& ss) {
  ss << "\"time_domain\" : [\n";

  for (std::size_t i = 0; i < elements_.size(); i++)
    if (i == elements_.size() - 1)
      ss << elements_[i] << "\n";
    else
      ss << elements_[i] << ",\n";

  ss << "]\n";
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_TIME_DOMAIN_HPP
