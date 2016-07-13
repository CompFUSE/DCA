// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef PHYS_LIBRARY_DOMAINS_TIME_AND_FREQUENCY_TIME_DOMAIN_H
#define PHYS_LIBRARY_DOMAINS_TIME_AND_FREQUENCY_TIME_DOMAIN_H

#include <cmath>
#include <string>
#include <vector>

#include "math_library/functional_transforms/domain_specifications/domain_specifications.hpp"

class time_domain {
public:
  const static int RULE = 1;
  const static int DIMENSION = 1;

  typedef double scalar_type;
  typedef double element_type;

  typedef math_algorithms::interval_dmn_1D_type dmn_specifications_type;
  
public:
  // const static int DIMENSION = 1;
  typedef time_domain parameter_type;  // --> used in the interpolation!

  static int& get_size();

  static std::string get_name();

  static std::vector<double>& get_elements();

  static double get_beta();

  static double get_volume();

  template <class stream_type>
  static void to_JSN(stream_type& ss);

  template <typename Reader>
  static void read(Reader& reader);

  template <typename Writer>
  static void write(Writer& writer);

  template <typename parameters_t>
  static void initialize(parameters_t& parameters);

  static void initialize_integration_domain(int level, std::vector<scalar_type>& weights,
                                            std::vector<element_type>& elements);

private:
  static int time_slices;
  static double beta;
};

int time_domain::time_slices = -1;
double time_domain::beta = 0;

int& time_domain::get_size() {
  static int size = -1;
  return size;
}

std::string time_domain::get_name() {
  static std::string name = "time-domain";
  return name;
}

std::vector<double>& time_domain::get_elements() {
  static std::vector<double> elements;
  return elements;
}

double time_domain::get_beta() {
  return (get_elements().back() + 1.e-10);
}

double time_domain::get_volume() {
  return (get_elements().back() + 1.e-10);
}

template <typename Writer>
void time_domain::write(Writer& writer) {
  writer.open_group(get_name());

  writer.execute("elements", get_elements());

  writer.close_group();
}

template <typename parameters_t>
void time_domain::initialize(parameters_t& parameters) {
  time_slices = parameters.get_sp_time_intervals();
  beta = parameters.get_beta();

  get_size() = 2 * (parameters.get_sp_time_intervals() + 1);

  get_elements().resize(get_size());

  for (int i = 0; i < get_size() / 2; i++) {
    get_elements()[i + get_size() / 2] = double(i) / (double(get_size()) / 2. - 1.) * beta;
    get_elements()[i] = -beta + double(i) / (double(get_size()) / 2. - 1.) * beta;
  }

  get_elements()[0] += 1.e-10;
  get_elements()[get_size() / 2 - 1] -= 1.e-10;
  get_elements()[get_size() / 2] += 1.e-10;
  get_elements()[get_size() / 2 - 1] -= 1.e-10;
}

template <class stream_type>
void time_domain::to_JSN(stream_type& ss) {
  ss << "\"time_domain\" : [\n";

  for (int i = 0; i < get_size(); i++)
    if (i == get_size() - 1)
      ss << get_elements()[i] << "\n";
    else
      ss << get_elements()[i] << ",\n";

  ss << "]\n";
}

void time_domain::initialize_integration_domain(int level, std::vector<scalar_type>& weights,
                                                std::vector<element_type>& elements) {
  int N = std::pow(2., level);

  weights.resize(0);
  elements.resize(0);

  for (int i = 0; i < get_size() / 2 - 1; i++) {
    element_type min = get_elements()[i];
    element_type max = get_elements()[i + 1];

    double delta = double(max - min) / double(N);

    for (int j = 0; j < N; j++) {
      elements.push_back(min + j * delta);

      weights.push_back(get_volume());
    }
  }

  for (int i = get_size() / 2; i < get_size() - 1; i++) {
    element_type min = get_elements()[i];
    element_type max = get_elements()[i + 1];

    double delta = double(max - min) / double(N);

    for (int j = 0; j < N; j++) {
      elements.push_back(min + j * delta);

      weights.push_back(get_volume());
    }
  }
}

#endif  // PHYS_LIBRARY_DOMAINS_TIME_AND_FREQUENCY_TIME_DOMAIN_H
