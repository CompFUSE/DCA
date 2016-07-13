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

#ifndef PHYS_LIBRARY_DOMAINS_QUANTUM_DOMAIN_POINT_GROUP_OPERATION_DMN_H
#define PHYS_LIBRARY_DOMAINS_QUANTUM_DOMAIN_POINT_GROUP_OPERATION_DMN_H

#include <cassert>
#include <complex>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

enum symmetry_group_level { UNIT_CELL, SUPER_CELL };
typedef symmetry_group_level symmetry_group_level_type;

std::string to_str(symmetry_group_level name) {
  switch (name) {
    case UNIT_CELL:
      return "UNIT_CELL";
      break;

    case SUPER_CELL:
      return "SUPER_CELL";
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

class point_group_symmetry_element {
public:
  /*!
   * should not be implemented. Just put here so we get a
   * linking error if the compiler wants to use it
   */
  point_group_symmetry_element();

  point_group_symmetry_element(int d);
  point_group_symmetry_element(const point_group_symmetry_element& copy_from_me);

  ~point_group_symmetry_element();

  void set_permutation(std::vector<int> p_vec);
  std::vector<int> get_permutation();

  void linear_transform(double* t0, double* t1);

  void transform(double* t0, double* t1);

  template <class stream_type>
  void to_JSON(stream_type& ss);

public:
  int DIMENSION;

  int ORDER;  // the order is defined such that O^N=1
  std::complex<double> PHASE;

  std::vector<int> P;

  double* O;
  double* t;
};

point_group_symmetry_element::point_group_symmetry_element(int d)
    : DIMENSION(d),

      ORDER(1),
      PHASE(1.),

      O(NULL),
      t(NULL) {
  O = new double[DIMENSION * DIMENSION];
  t = new double[DIMENSION];

  memset(t, 0, DIMENSION * sizeof(double));
  memset(O, 0, DIMENSION * DIMENSION * sizeof(double));

  for (int j = 0; j < DIMENSION; ++j)
    for (int i = 0; i < DIMENSION; ++i)
      O[i + j * DIMENSION] = i == j ? 1 : 0;
}

point_group_symmetry_element::point_group_symmetry_element(const point_group_symmetry_element& other)
    : DIMENSION(other.DIMENSION),

      ORDER(other.ORDER),
      PHASE(other.PHASE),

      O(NULL),
      t(NULL) {
  P = other.P;

  O = new double[DIMENSION * DIMENSION];
  t = new double[DIMENSION];

  memcpy(t, other.t, DIMENSION * sizeof(double));
  memcpy(O, other.O, DIMENSION * DIMENSION * sizeof(double));
}

point_group_symmetry_element::~point_group_symmetry_element() {
  delete[] O;
  delete[] t;
}

void point_group_symmetry_element::set_permutation(std::vector<int> p_vec) {
  P = p_vec;
}

std::vector<int> point_group_symmetry_element::get_permutation() {
  return P;
}

void point_group_symmetry_element::linear_transform(double* t0, double* t1) {
  for (int i = 0; i < DIMENSION; ++i)
    t1[i] = 0;

  for (int j = 0; j < DIMENSION; ++j)
    for (int i = 0; i < DIMENSION; ++i)
      t1[i] += O[i + DIMENSION * j] * t0[j];
}

void point_group_symmetry_element::transform(double* t0, double* t1) {
  for (int i = 0; i < DIMENSION; ++i)
    t1[i] = 0;

  for (int j = 0; j < DIMENSION; ++j)
    for (int i = 0; i < DIMENSION; ++i)
      t1[i] += O[i + DIMENSION * j] * t0[j];

  for (int i = 0; i < DIMENSION; ++i)
    t1[i] += t[i];
}

template <class stream_type>
void point_group_symmetry_element::to_JSON(stream_type& ss) {
  ss << "\"O\" : [\n";

  for (int i = 0; i < DIMENSION; ++i) {
    ss << "[";
    for (int j = 0; j < DIMENSION; ++j) {
      if (j == DIMENSION - 1)
        ss << O[i + j * DIMENSION] << "]";
      else
        ss << O[i + j * DIMENSION] << ", ";
    }

    if (i == DIMENSION - 1)
      ss << "\n],\n";
    else
      ss << ",\n";
  }

  ss << "\"t\" : [ ";
  for (int i = 0; i < DIMENSION; ++i) {
    if (i == DIMENSION - 1)
      ss << t[i] << "]";
    else
      ss << t[i] << ", ";
  }
  ss << ",\n";

  ss << "\"P\" : [ ";
  if (P.size() == 0)
    ss << "-1]";
  else
    for (size_t i = 0; i < P.size(); ++i) {
      if (i == P.size() - 1)
        ss << P[i] << "]";
      else
        ss << P[i] << ", ";
    }

  ss << "\n\n\n";
}

template <symmetry_group_level_type symmetry_group_level, typename base_cluster_type>
class point_group_symmetry_domain {
public:
  typedef point_group_symmetry_element element_type;

public:
  static int DIMENSION;

public:
  static int& get_size();
  static std::string get_name();
  static std::vector<element_type>& get_elements();

  static bool is_initialized();

  static std::vector<double> transform(int i, std::vector<double> vec);

  template <typename Writer>
  static void write(Writer& writer);

  template <class stream_type>
  static void to_JSON(stream_type& ss);
};

template <symmetry_group_level_type symmetry_group_level, typename base_cluster_type>
int point_group_symmetry_domain<symmetry_group_level, base_cluster_type>::DIMENSION = -1;

template <symmetry_group_level_type symmetry_group_level, typename base_cluster_type>
int& point_group_symmetry_domain<symmetry_group_level, base_cluster_type>::get_size() {
  static int size = 0;
  return size;
}

template <symmetry_group_level_type symmetry_group_level, typename base_cluster_type>
std::string point_group_symmetry_domain<symmetry_group_level, base_cluster_type>::get_name() {
  static std::string name = "point-group-symmetry-domain (" + to_str(symmetry_group_level) + ", " +
                            base_cluster_type::get_name() + ")";

  return name;
}

template <symmetry_group_level_type symmetry_group_level, typename base_cluster_type>
std::vector<point_group_symmetry_element>& point_group_symmetry_domain<
    symmetry_group_level, base_cluster_type>::get_elements() {
  static std::vector<point_group_symmetry_element> v(get_size(),
                                                     point_group_symmetry_element(DIMENSION));
  return v;
}

template <symmetry_group_level_type symmetry_group_level, typename base_cluster_type>
bool point_group_symmetry_domain<symmetry_group_level, base_cluster_type>::is_initialized() {
  static bool initialized = false;
  return initialized;
}

template <symmetry_group_level_type symmetry_group_level, typename base_cluster_type>
std::vector<double> point_group_symmetry_domain<symmetry_group_level, base_cluster_type>::transform(
    int i, std::vector<double> vec) {
  assert(i > -1);
  assert(i < get_size());

  std::vector<double> tmp(DIMENSION, 0);

  get_elements()[i].transform(&vec[0], &tmp[0]);

  return tmp;
}

template <symmetry_group_level_type symmetry_group_level, typename base_cluster_type>
template <typename Writer>
void point_group_symmetry_domain<symmetry_group_level, base_cluster_type>::write(Writer& /*writer*/) {
  //   std::vector<double>               T(DIMENSION);
  //   std::vector<std::vector<double> > O(DIMENSION, std::vector<double>(DIMENSION));

  //   writer.open_group(get_name());

  //   writer.close_group();
}

template <symmetry_group_level_type symmetry_group_level, typename base_cluster_type>
template <class stream_type>
void point_group_symmetry_domain<symmetry_group_level, base_cluster_type>::to_JSON(stream_type& ss) {
  ss << "\"pointgroup symmetry domain ";

  if (symmetry_group_level == UNIT_CELL)
    ss << "unit-cell";

  if (symmetry_group_level == SUPER_CELL)
    ss << "super-cell";

  ss << "\" : {\n";

  for (int l = 0; l < get_size(); ++l) {
    ss << "\"" << l << "\" : \n{\n";

    get_elements()[l].to_JSON(ss);

    if (l == get_size() - 1)
      ss << "\n}\n";
    else
      ss << "\n},\n";
  }
  ss << "}\n";
}

#endif  // PHYS_LIBRARY_DOMAINS_QUANTUM_DOMAIN_POINT_GROUP_OPERATION_DMN_H
