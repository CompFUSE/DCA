// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// The domain contained in this file describes the frequency exchange inside a two particle Green's
// Function.

#ifndef DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_FREQUENCY_EXCHANGE_DOMAIN_HPP
#define DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_FREQUENCY_EXCHANGE_DOMAIN_HPP

#include <cassert>
#include <string>
#include <vector>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

class FrequencyExchangeDomain {
public:
  using scalar_type = int;
  using element_type = int;

  // Initializes size and elements.
  template <class Parameters>
  static void initialize(const Parameters& parameters);

  // Returns the number of computed frequency exchanges.
  // Precondition: the domain is initialized.
  static int get_size() {
    assert(initialized_);
    return elements_.size();
  }

  // Contains all indices relative to the exchanged frequencies inside the two particle
  // computation.
  // Precondition: the domain is initialized.
  static const std::vector<int>& get_elements() {
    assert(initialized_);
    return elements_;
  }

  // Returns the number of additional frequencies to store in G1.
  static int extensionSize();

  static bool isInitialized() {
    return initialized_;
  }

  static const std::string& get_name() {
    const static std::string name = "Frequency exchange domain.";
    return name;
  }

  template <class Writer>
  static void write(Writer& writer);

private:
  static std::vector<int> elements_;
  static bool initialized_;
};
std::vector<int> FrequencyExchangeDomain::elements_;
bool FrequencyExchangeDomain::initialized_ = false;

template <class Parameters>
void FrequencyExchangeDomain::initialize(const Parameters& parameters) {
  if (parameters.compute_all_transfers()) {
    assert(parameters.get_four_point_frequency_transfer() > 0);
    elements_.resize(parameters.get_four_point_frequency_transfer() + 1);
    int idx_value = 0;
    for (int& elem : elements_)
      elem = idx_value++;
  }

  else {
    elements_ = std::vector<int>{parameters.get_four_point_frequency_transfer()};
  }

  initialized_ = true;
}

int FrequencyExchangeDomain::extensionSize() {
  auto compute_extension = [] {
    if (!initialized_)
      throw(std::logic_error("The frequency exchange domain was not initialized."));
    int size = 0;
    for (auto el : elements_)
      size = std::max(size, std::abs(el));
    return size;
  };
  static int extension_size = compute_extension();
  return extension_size;
}

template <class Writer>
void FrequencyExchangeDomain::write(Writer& writer) {
  writer.open_group(get_name());
  writer.execute("element_indices", elements_);
  writer.close_group();
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_FREQUENCY_EXCHANGE_DOMAIN_HPP
