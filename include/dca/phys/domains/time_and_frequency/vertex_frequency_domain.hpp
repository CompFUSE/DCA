// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements various types and orderings of the matsubara frequency domain for the
// vertex function via templates.

#ifndef DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_VERTEX_FREQUENCY_DOMAIN_HPP
#define DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_VERTEX_FREQUENCY_DOMAIN_HPP

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_exchange_domain.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_frequency_name.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <VERTEX_FREQUENCY_NAME NAME>
class vertex_frequency_domain {
public:
  const static int DIMENSION = 1;

  typedef double scalar_type;
  typedef double element_type;

  static int& get_size() {
    static int size;
    return size;
  }

  static std::string get_name() {
    static std::string name = "vertex-frequency-domain (" + to_str(NAME) + ")";
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
    static std::vector<element_type> elements;
    return elements;
  }

  static std::vector<int>& get_corresponding_frequency_domain_index() {
    static std::vector<int> elements;
    return elements;
  }

  template <typename Writer>
  static void write(Writer& writer);

  template <typename parameters_t>
  static void initialize(parameters_t& parameters);
};

template <VERTEX_FREQUENCY_NAME NAME>
template <typename Writer>
void vertex_frequency_domain<NAME>::write(Writer& writer) {
  writer.open_group(get_name());
  writer.execute("elements", get_elements());
  writer.close_group();
}

template <>
template <typename parameters_t>
void vertex_frequency_domain<COMPACT>::initialize(parameters_t& parameters) {
  get_basis()[0] = (2. * M_PI) / parameters.get_beta();
  get_inverse_basis()[0] = parameters.get_beta() / (2. * M_PI);

  get_size() = 2 * parameters.get_four_point_fermionic_frequencies();  // wn_c();

  get_elements().resize(get_size());

  for (int l = 0; l < parameters.get_four_point_fermionic_frequencies(); l++) {
    get_elements()[get_size() / 2 + l] = M_PI / parameters.get_beta() * (1 + 2 * l);
    get_elements()[get_size() / 2 - 1 - l] = -M_PI / parameters.get_beta() * (1 + 2 * l);
  }

  get_corresponding_frequency_domain_index().resize(get_size(), -1);

  const std::vector<double>& wn = frequency_domain::get_elements();

  for (int i = 0; i < get_size(); i++)
    for (size_t j = 0; j < wn.size(); j++)
      if (std::fabs(wn[j] - get_elements()[i]) < 1.e-6)
        get_corresponding_frequency_domain_index()[i] = j;

  for (int i = 0; i < get_size(); i++)
    if (get_corresponding_frequency_domain_index()[i] == -1 ||
        std::fabs(wn[get_corresponding_frequency_domain_index()[i]] - get_elements()[i]) > 1.e-6)
      throw std::logic_error(__FUNCTION__);
}

template <>
template <typename parameters_t>
void vertex_frequency_domain<COMPACT_POSITIVE>::initialize(parameters_t& parameters) {
  get_basis()[0] = (2. * M_PI) / parameters.get_beta();
  get_inverse_basis()[0] = parameters.get_beta() / (2. * M_PI);

  get_size() = parameters.get_four_point_fermionic_frequencies();  // wn_c();

  get_elements().resize(get_size());

  for (int l = 0; l < parameters.get_four_point_fermionic_frequencies(); l++)
    get_elements()[l] = M_PI / parameters.get_beta() * (1 + 2 * l);

  get_corresponding_frequency_domain_index().resize(get_size(), -1);

  const std::vector<double>& wn = frequency_domain::get_elements();

  for (int i = 0; i < get_size(); i++)
    for (size_t j = 0; j < wn.size(); j++)
      if (std::fabs(wn[j] - get_elements()[i]) < 1.e-6)
        get_corresponding_frequency_domain_index()[i] = j;

  for (int i = 0; i < get_size(); i++)
    if (get_corresponding_frequency_domain_index()[i] == -1 ||
        std::fabs(wn[get_corresponding_frequency_domain_index()[i]] - get_elements()[i]) > 1.e-6)
      throw std::logic_error(__FUNCTION__);

  assert(get_elements().back() == vertex_frequency_domain<COMPACT>::get_elements().back());
}

template <>
template <typename parameters_t>
void vertex_frequency_domain<EXTENDED>::initialize(parameters_t& parameters) {
  get_basis()[0] = (2. * M_PI) / parameters.get_beta();
  get_inverse_basis()[0] = parameters.get_beta() / (2. * M_PI);

  if(!FrequencyExchangeDomain::isInitialized())
    FrequencyExchangeDomain::initialize(parameters);

  get_size() = 2 * (parameters.get_four_point_fermionic_frequencies() +
                    FrequencyExchangeDomain::extensionSize());

  get_elements().resize(get_size());

  for (int l = 0; l < get_size() / 2; l++) {
    get_elements()[get_size() / 2 + l] = M_PI / parameters.get_beta() * (1 + 2 * l);
    get_elements()[get_size() / 2 - 1 - l] = -M_PI / parameters.get_beta() * (1 + 2 * l);
  }

  get_corresponding_frequency_domain_index().resize(get_size(), -1);

  const std::vector<double>& wn = frequency_domain::get_elements();

  for (int i = 0; i < get_size(); i++)
    for (size_t j = 0; j < wn.size(); j++)
      if (std::fabs(wn[j] - get_elements()[i]) < 1.e-6)
        get_corresponding_frequency_domain_index()[i] = j;

  for (int i = 0; i < get_size(); i++)
    if (get_corresponding_frequency_domain_index()[i] == -1 ||
        std::fabs(wn[get_corresponding_frequency_domain_index()[i]] - get_elements()[i]) > 1.e-6)
      throw std::logic_error(__FUNCTION__);
}

template <>
template <typename parameters_t>
void vertex_frequency_domain<EXTENDED_POSITIVE>::initialize(parameters_t& parameters) {
  get_basis()[0] = (2. * M_PI) / parameters.get_beta();
  get_inverse_basis()[0] = parameters.get_beta() / (2. * M_PI);

  if(!FrequencyExchangeDomain::isInitialized())
    FrequencyExchangeDomain::initialize(parameters);

  get_size() =
      parameters.get_four_point_fermionic_frequencies() + FrequencyExchangeDomain::extensionSize();

  get_elements().resize(get_size());

  for (int l = 0; l < get_size(); l++)
    get_elements()[l] = M_PI / parameters.get_beta() * (1 + 2 * l);

  get_corresponding_frequency_domain_index().resize(get_size(), -1);

  const std::vector<double>& wn = frequency_domain::get_elements();

  for (int i = 0; i < get_size(); i++)
    for (size_t j = 0; j < wn.size(); j++)
      if (std::fabs(wn[j] - get_elements()[i]) < 1.e-6)
        get_corresponding_frequency_domain_index()[i] = j;

  for (int i = 0; i < get_size(); i++)
    if (get_corresponding_frequency_domain_index()[i] == -1 ||
        std::fabs(wn[get_corresponding_frequency_domain_index()[i]] - get_elements()[i]) > 1.e-6)
      throw std::logic_error(__FUNCTION__);

  assert(get_elements().back() == vertex_frequency_domain<EXTENDED>::get_elements().back());
}

template <>
template <typename parameters_t>
void vertex_frequency_domain<EXTENDED_BOSONIC>::initialize(parameters_t& parameters) {
  get_basis()[0] = (2. * M_PI) / parameters.get_beta();
  get_inverse_basis()[0] = parameters.get_beta() / (2. * M_PI);

  get_size() = 2 * parameters.get_hts_bosonic_frequencies() + 1;

  get_elements().resize(get_size(), -2. * M_PI / parameters.get_beta() * int(get_size() / 2));

  for (int l = 0; l < get_size(); l++) {
    get_elements()[l] += l * 2. * M_PI / parameters.get_beta();
    // cout << get_elements()[l] << endl;
  }
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_VERTEX_FREQUENCY_DOMAIN_HPP
