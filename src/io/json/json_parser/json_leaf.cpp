// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements json_leaf.hpp.

#include "dca/io/json/json_parser/json_leaf.hpp"
#include <cstring>

namespace dca {
namespace io {
// dca::io::

JSON_leaf::JSON_leaf()
    : N(0),
      LD(0),

      bool_ptr(NULL),
      char_ptr(NULL),
      int_ptr(NULL),
      double_ptr(NULL) {}

JSON_leaf::JSON_leaf(int n)
    : N(n),
      LD(n),

      bool_ptr(NULL),
      char_ptr(NULL),
      int_ptr(NULL),
      double_ptr(NULL) {
  bool_ptr = new bool[LD];
  char_ptr = new char[LD];
}

JSON_leaf::JSON_leaf(int n, int ld)
    : N(n),
      LD(ld),

      bool_ptr(NULL),
      char_ptr(NULL),
      int_ptr(NULL),
      double_ptr(NULL) {}

JSON_leaf::~JSON_leaf() {
  if (bool_ptr != NULL)
    delete[] bool_ptr;

  if (int_ptr != NULL)
    delete[] int_ptr;

  if (char_ptr != NULL)
    delete[] char_ptr;

  if (double_ptr != NULL)
    delete[] double_ptr;
}

template <typename ss>
void JSON_leaf::print(ss& ss_obj) {
  for (size_t i = 0; i < N; i++)
    ss_obj << char_ptr[i];
}

JSON_leaf& JSON_leaf::operator=(bool rhs) {
  if (bool_ptr == NULL) {
    N = 1;
    bool_ptr = new bool[N];
  }

  bool_ptr[0] = rhs;
  return *this;
}

JSON_leaf& JSON_leaf::operator=(char rhs) {
  if (char_ptr == NULL) {
    N = 1;
    char_ptr = new char[N];
  }

  char_ptr[0] = rhs;
  return *this;
}

JSON_leaf& JSON_leaf::operator=(std::string rhs) {
  N = rhs.size();

  if (char_ptr != NULL)
    delete[] char_ptr;

  char_ptr = new char[N];

  std::memcpy(char_ptr, &rhs[0], N * sizeof(char));

  return *this;
}

JSON_leaf& JSON_leaf::operator=(int rhs) {
  if (int_ptr == NULL) {
    N = 1;
    int_ptr = new int[N];
  }

  int_ptr[0] = rhs;
  return *this;
}

JSON_leaf& JSON_leaf::operator=(double rhs) {
  if (double_ptr == NULL) {
    N = 1;
    double_ptr = new double[N];
  }

  double_ptr[0] = rhs;
  return *this;
}

}  // io
}  // dca
