// Copyright (C) 2026 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Anirudha Mirmira
//
// Runtime-sized index domain enumerating the disorder configurations of a DCA step. It is
// used to assemble the combined "disorder-randoms" function written when the disorder
// "dump-randoms" parameter is set. Its size is established by initialize(num_configurations)
// and is not part of the global domain initialization.

#ifndef DCA_PHYS_DOMAINS_CLUSTER_DISORDER_CONFIGURATION_DOMAIN_HPP
#define DCA_PHYS_DOMAINS_CLUSTER_DISORDER_CONFIGURATION_DOMAIN_HPP

#include <numeric>
#include <string>
#include <vector>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

class DisorderConfigurationDomain {
public:
  using scalar_type = int;
  using element_type = int;

  // Set the number of disorder configurations. Call before the domain is used.
  static void initialize(int num_configurations);

  // Returns the number of disorder configurations.
  static int get_size() {
    return elements_.size();
  }

  // The configuration indices 0 .. get_size()-1.
  static const std::vector<int>& get_elements() {
    return elements_;
  }

  static const std::string& get_name() {
    const static std::string name = "Disorder configuration domain.";
    return name;
  }

  template <class Writer>
  static void write(Writer& writer);

private:
  inline static std::vector<int> elements_{};
};

inline void DisorderConfigurationDomain::initialize(int num_configurations) {
  elements_.resize(num_configurations);
  std::iota(elements_.begin(), elements_.end(), 0);
}

template <class Writer>
void DisorderConfigurationDomain::write(Writer& writer) {
  writer.open_group(get_name());
  writer.execute("element_indices_", elements_);
  writer.close_group();
}

}  // namespace domains
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DOMAINS_CLUSTER_DISORDER_CONFIGURATION_DOMAIN_HPP
