// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class implements various types and orderings of the matsubara-time domain via templates.

#ifndef PHYS_LIBRARY_DOMAINS_TIME_AND_FREQUENCY_VERTEX_TIME_DOMAIN_H
#define PHYS_LIBRARY_DOMAINS_TIME_AND_FREQUENCY_VERTEX_TIME_DOMAIN_H

#include <stdexcept>
#include <string>
#include <vector>

namespace DCA {
enum VERTEX_TIME_NAME {
  SP_TIME_DOMAIN,
  SP_TIME_DOMAIN_POSITIVE,
  SP_TIME_DOMAIN_LEFT_ORIENTED,
  SP_TIME_DOMAIN_LEFT_ORIENTED_POSITIVE,
  TP_TIME_DOMAIN,
  TP_TIME_DOMAIN_POSITIVE
};

std::string to_str(VERTEX_TIME_NAME NAME) {
  switch (NAME) {
    case SP_TIME_DOMAIN:
      return "SP_TIME_DOMAIN";
      break;

    case SP_TIME_DOMAIN_POSITIVE:
      return "SP_TIME_DOMAIN_POSITIVE";
      break;

    case SP_TIME_DOMAIN_LEFT_ORIENTED:
      return "SP_TIME_DOMAIN_LEFT_ORIENTED";
      break;

    case SP_TIME_DOMAIN_LEFT_ORIENTED_POSITIVE:
      return "SP_TIME_DOMAIN_LEFT_ORIENTED_POSITIVE";
      break;

    case TP_TIME_DOMAIN:
      return "TP_TIME_DOMAIN";
      break;

    case TP_TIME_DOMAIN_POSITIVE:
      return "TP_TIME_DOMAIN_POSITIVE";
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }

  return "NONE";
}

template <VERTEX_TIME_NAME NAME>
class vertex_time_domain {
public:
  const static int RULE = 1;
  const static int DIMENSION = 1;

  typedef double scalar_type;
  typedef double element_type;

  typedef vertex_time_domain<NAME> parameter_type;

public:
  static int& get_size();
  static std::string get_name();

  static scalar_type* get_basis();
  static scalar_type* get_inverse_basis();

  static scalar_type* get_super_basis();
  static scalar_type* get_inverse_super_basis();

  static std::vector<element_type>& get_elements();

  template <typename Reader>
  static void read(Reader& reader);

  template <typename Writer>
  static void write(Writer& writer);

  template <typename parameters_t>
  static void initialize(parameters_t& parameters);
};

template <VERTEX_TIME_NAME NAME>
int& vertex_time_domain<NAME>::get_size() {
  static int size = -1;
  return size;
}

template <VERTEX_TIME_NAME NAME>
std::string vertex_time_domain<NAME>::get_name() {
  static std::string name = to_str(NAME);
  return name;
}

template <VERTEX_TIME_NAME NAME>
typename vertex_time_domain<NAME>::scalar_type* vertex_time_domain<NAME>::get_basis() {
  static scalar_type basis[DIMENSION];
  return basis;
}

template <VERTEX_TIME_NAME NAME>
typename vertex_time_domain<NAME>::scalar_type* vertex_time_domain<NAME>::get_inverse_basis() {
  static scalar_type inverse_basis[DIMENSION];
  return inverse_basis;
}

template <VERTEX_TIME_NAME NAME>
typename vertex_time_domain<NAME>::scalar_type* vertex_time_domain<NAME>::get_super_basis() {
  static scalar_type super_basis[DIMENSION];
  return super_basis;
}

template <VERTEX_TIME_NAME NAME>
typename vertex_time_domain<NAME>::scalar_type* vertex_time_domain<NAME>::get_inverse_super_basis() {
  static scalar_type inverse_super_basis[DIMENSION];
  return inverse_super_basis;
}

template <VERTEX_TIME_NAME NAME>
std::vector<typename vertex_time_domain<NAME>::element_type>& vertex_time_domain<NAME>::get_elements() {
  static std::vector<element_type> elements;
  return elements;
}

template <VERTEX_TIME_NAME NAME>
template <typename Writer>
void vertex_time_domain<NAME>::write(Writer& writer) {
  writer.open_group(get_name());

  writer.execute("elements", get_elements());

  writer.close_group();
}

template <VERTEX_TIME_NAME NAME>
template <typename parameters_t>
void vertex_time_domain<NAME>::initialize(parameters_t& /*parameters*/) {}

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
  get_elements()[get_size() / 2 - 1] -= 1.e-10;
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
  int time_slices = parameters.get_tp_time_intervals();
  double beta = parameters.get_beta();

  get_size() = 2 * (parameters.get_tp_time_intervals() + 1);

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
  get_elements()[get_size() / 2 - 1] -= 1.e-10;
}

template <>
template <typename parameters_t>
void vertex_time_domain<TP_TIME_DOMAIN_POSITIVE>::initialize(parameters_t& parameters) {
  int time_slices = parameters.get_tp_time_intervals();
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
}

#endif  // PHYS_LIBRARY_DOMAINS_TIME_AND_FREQUENCY_VERTEX_TIME_DOMAIN_H
