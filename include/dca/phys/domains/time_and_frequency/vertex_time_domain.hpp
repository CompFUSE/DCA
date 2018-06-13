// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements various types and orderings of the matsubara time domain for the vertex
// function via templates.

#ifndef DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_VERTEX_TIME_DOMAIN_HPP
#define DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_VERTEX_TIME_DOMAIN_HPP

#include <string>
#include <vector>

#include "dca/phys/domains/time_and_frequency/vertex_time_name.hpp"

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

template <VERTEX_TIME_NAME NAME>
class vertex_time_domain {
public:
  const static int RULE = 1;
  const static int DIMENSION = 1;

  typedef double scalar_type;
  typedef double element_type;

  typedef vertex_time_domain<NAME> parameter_type;

  static int& get_size() {
    static int size = -1;
    return size;
  }

  static std::string get_name() {
    static std::string name = to_str(NAME);
    return name;
  }

  static scalar_type* get_basis() {
    static scalar_type basis[DIMENSION];
    return basis;
  }

  static scalar_type* get_inverse_basis() {
    static scalar_type inverse_basis[DIMENSION];
    return inverse_basis;
  }

  static scalar_type* get_super_basis() {
    static scalar_type super_basis[DIMENSION];
    return super_basis;
  }

  static scalar_type* get_inverse_super_basis() {
    static scalar_type inverse_super_basis[DIMENSION];
    return inverse_super_basis;
  }

  static std::vector<element_type>& get_elements() {
    static std::vector<element_type> elements;
    return elements;
  }

  template <typename Writer>
  static void write(Writer& writer);

  template <typename parameters_t>
  static void initialize(parameters_t& parameters);
};

template <VERTEX_TIME_NAME NAME>
template <typename Writer>
void vertex_time_domain<NAME>::write(Writer& writer) {
  writer.open_group(get_name());
  writer.execute("elements", get_elements());
  writer.close_group();
}

template <>
template <typename parameters_t>
void vertex_time_domain<SP_TIME_DOMAIN>::initialize(parameters_t& parameters) {
  int time_slices = parameters.get_sp_time_intervals();
  double beta = parameters.get_beta();

  get_size() = 2 * (parameters.get_sp_time_intervals() + 1);

  get_basis()[0] = beta / double(time_slices);
  get_inverse_basis()[0] = double(time_slices) / beta;

  get_super_basis()[0] = 2 * beta;
  get_inverse_super_basis()[0] = 1. / (2 * beta);

  get_elements().resize(get_size());

  for (int i = 0; i < get_size() / 2; i++) {
    get_elements()[i + get_size() / 2] = double(i) / double(time_slices) * beta;
    get_elements()[i] = -beta + double(i) / double(time_slices) * beta;
  }

  get_elements()[0] += 1.e-10;
  get_elements()[get_size() / 2 - 1] -= 1.e-10;
  get_elements()[get_size() / 2] += 1.e-10;
  get_elements()[get_size() / 2 - 1] -= 1.e-10;  // TODO: Typo? (see time_domain.hpp)
}

template <>
template <typename parameters_t>
void vertex_time_domain<SP_TIME_DOMAIN_POSITIVE>::initialize(parameters_t& parameters) {
  int time_slices = parameters.get_sp_time_intervals();
  double beta = parameters.get_beta();

  get_size() = time_slices + 1;

  get_basis()[0] = beta / double(time_slices);
  get_inverse_basis()[0] = double(time_slices) / beta;

  get_super_basis()[0] = beta;
  get_inverse_super_basis()[0] = 1. / beta;

  get_elements().resize(get_size());

  for (int i = 0; i < get_size(); i++)
    get_elements()[i] = double(i) / double(time_slices) * beta;

  get_elements()[0] += 1.e-10;
  get_elements()[time_slices - 1] -= 1.e-10;
}

template <>
template <typename parameters_t>
void vertex_time_domain<SP_TIME_DOMAIN_LEFT_ORIENTED>::initialize(parameters_t& parameters) {
  int time_slices = parameters.get_sp_time_intervals();
  double beta = parameters.get_beta();

  get_size() = 2 * time_slices;

  get_basis()[0] = beta / double(time_slices);
  get_inverse_basis()[0] = double(time_slices) / beta;

  get_super_basis()[0] = 2 * beta;
  get_inverse_super_basis()[0] = 1. / (2 * beta);

  get_elements().resize(get_size());

  for (int i = 0; i < get_size() / 2; i++) {
    get_elements()[i + get_size() / 2] = double(i) / double(time_slices) * beta;
    get_elements()[i] = -beta + double(i) / double(time_slices) * beta;
  }

  get_elements()[0] += 1.e-10;
  get_elements()[get_size() / 2] += 1.e-10;
}

template <>
template <typename parameters_t>
void vertex_time_domain<SP_TIME_DOMAIN_LEFT_ORIENTED_POSITIVE>::initialize(parameters_t& parameters) {
  int time_slices = parameters.get_sp_time_intervals();
  double beta = parameters.get_beta();

  get_size() = time_slices;

  get_basis()[0] = beta / double(time_slices);
  get_inverse_basis()[0] = double(time_slices) / beta;

  get_super_basis()[0] = beta;
  get_inverse_super_basis()[0] = 1. / beta;

  get_elements().resize(get_size());

  for (int i = 0; i < get_size(); i++)
    get_elements()[i] = double(i) / double(time_slices) * beta;

  get_elements()[0] += 1.e-10;
}

template <>
template <typename parameters_t>
void vertex_time_domain<TP_TIME_DOMAIN>::initialize(parameters_t& parameters) {
  int time_slices = parameters.get_time_intervals_for_time_measurements();
  double beta = parameters.get_beta();

  get_size() = 2 * (parameters.get_time_intervals_for_time_measurements() + 1);

  get_basis()[0] = beta / double(time_slices);
  get_inverse_basis()[0] = double(time_slices) / beta;

  get_super_basis()[0] = 2 * beta;
  get_inverse_super_basis()[0] = 1. / (2 * beta);

  get_elements().resize(get_size());

  for (int i = 0; i < get_size() / 2; i++) {
    get_elements()[i + get_size() / 2] = double(i) / double(time_slices) * beta;
    get_elements()[i] = -beta + double(i) / double(time_slices) * beta;
  }

  get_elements()[0] += 1.e-10;
  get_elements()[get_size() / 2 - 1] -= 1.e-10;
  get_elements()[get_size() / 2] += 1.e-10;
  get_elements()[get_size() / 2 - 1] -= 1.e-10;  // TODO: Typo?
}

template <>
template <typename parameters_t>
void vertex_time_domain<TP_TIME_DOMAIN_POSITIVE>::initialize(parameters_t& parameters) {
  int time_slices = parameters.get_time_intervals_for_time_measurements();
  double beta = parameters.get_beta();

  get_size() = time_slices + 1;

  get_basis()[0] = beta / double(time_slices);
  get_inverse_basis()[0] = double(time_slices) / beta;

  get_super_basis()[0] = beta;
  get_inverse_super_basis()[0] = 1. / beta;

  get_elements().resize(get_size());

  for (int i = 0; i < get_size(); i++) {
    get_elements()[i] = double(i) / double(time_slices) * beta;
  }

  get_elements()[0] += 1.e-10;
  get_elements()[get_size() - 1] -= 1.e-10;

  /*
  cout.precision(16);
  cout << get_name() << "\n";
  for(int i=0; i<get_size(); i++){
    cout << i << "\t" << get_elements()[i] << "\n";
  }
  cout << "\n";
  */
}

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_TIME_AND_FREQUENCY_VERTEX_TIME_DOMAIN_HPP
