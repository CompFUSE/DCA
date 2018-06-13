// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides an element of a point group symmetry.

#ifndef DCA_PHYS_DOMAINS_QUANTUM_POINT_GROUP_SYMMETRY_ELEMENT_HPP
#define DCA_PHYS_DOMAINS_QUANTUM_POINT_GROUP_SYMMETRY_ELEMENT_HPP

#include <complex>
#include <vector>

namespace dca {
namespace phys {
namespace domains {
// dca::phys::domains::

class point_group_symmetry_element {
public:
  point_group_symmetry_element(int d);
  point_group_symmetry_element(const point_group_symmetry_element& copy_from_me);

  ~point_group_symmetry_element();

  void set_permutation(std::vector<int> p_vec) {
    P = p_vec;
  }
  std::vector<int> get_permutation() {
    return P;
  }

  void linear_transform(double* t0, double* t1);

  void transform(double* t0, double* t1);

  template <class stream_type>
  void to_JSON(stream_type& ss);

  int DIMENSION;

  int ORDER;  // the order is defined such that O^N=1
  std::complex<double> PHASE;

  std::vector<int> P;

  double* O;
  double* t;
};

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

}  // domains
}  // phys
}  // dca

#endif  // DCA_PHYS_DOMAINS_QUANTUM_POINT_GROUP_SYMMETRY_ELEMENT_HPP
