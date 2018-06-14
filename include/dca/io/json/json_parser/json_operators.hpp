// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides operators for the JSON parser.

#ifndef DCA_IO_JSON_JSON_PARSER_JSON_OPERATORS_HPP
#define DCA_IO_JSON_JSON_PARSER_JSON_OPERATORS_HPP

#include <complex>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "dca/io/json/json_parser/parse_buffer.hpp"
#include "dca/io/json/json_parser/whatever.hpp"

namespace dca {
namespace io {
// dca::io::

bool& operator<=(bool& lhs, const Whatever& w);
int& operator<=(int& lhs, const Whatever& w);
float& operator<=(float& lhs, const Whatever& w);
double& operator<=(double& lhs, const Whatever& w);
std::string& operator<=(std::string& lhs, const Whatever& w);

template <typename T>
std::complex<T>& operator<=(std::complex<T>& lhs, const Whatever& w) {
  switch (w.type) {
    case WHATEVER_VECTOR: {
      T v0;
      v0 <= w.whateverVector[0];
      T v1;
      v1 <= w.whateverVector[1];
      std::complex<T> result(v0, v1);
      lhs = result;
      return lhs;
    }

    default: { throw std::logic_error(__FUNCTION__); }
  }
}

template <typename T>
std::vector<T>& operator<=(std::vector<T>& lhs, const Whatever& w) {
  switch (w.type) {
    case WHATEVER_VECTOR: {
      lhs.resize(w.whateverVector.size());
      for (size_t i = 0; i < lhs.size(); i++)
        lhs[i] <= w.whateverVector[i];
      return lhs;
    }

    default: { throw std::logic_error(__FUNCTION__); }
  }
}

template <typename T>
std::vector<std::vector<T>>& operator<=(std::vector<std::vector<T>>& lhs, const Whatever& w) {
  switch (w.type) {
    case WHATEVER_VECTOR: {
      lhs.resize(w.whateverVector.size(), std::vector<T>(w[0].size()));
      for (size_t i = 0; i < lhs.size(); i++)
        lhs[i] <= w.whateverVector[i];
      return lhs;
    }

    default: { throw std::logic_error(__FUNCTION__); }
  }
}

template <typename T>
std::map<std::string, T>& operator<=(std::map<std::string, T>& lhs, const Whatever& w) {
  switch (w.type) {
    case WHATEVER_MAP: {
      typedef std::map<std::wstring, Whatever> WhateverMapType;
      typedef typename WhateverMapType::const_iterator itr_type;

      lhs.clear();

      for (itr_type itr = w.whateverMap.begin(); itr != w.whateverMap.end(); itr++) {
        const std::pair<std::wstring, Whatever>& pair(*itr);
        std::string key(pair.first.begin(), pair.first.end());
        lhs[key] <= pair.second;
      }

      return lhs;
    }

    default: { throw std::logic_error(__FUNCTION__); }
  }
}

template <typename T>
ParseBuffer& operator>>(ParseBuffer& buffer, T& value) {
  std::wistringstream is(buffer.str());
  is >> value;
  return buffer;
}

}  // io
}  // dca

#endif  // DCA_IO_JSON_JSON_PARSER_JSON_OPERATORS_HPP
