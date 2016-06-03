// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// Description

#ifndef COMP_LIBRARY_IO_LIBRARY_JSON_JSON_PARSER_JSON_OPERATORS_H
#define COMP_LIBRARY_IO_LIBRARY_JSON_JSON_PARSER_JSON_OPERATORS_H

#include <complex>

#include "comp_library/IO_library/JSON/JSON_PARSER/parse_buffer.h"
#include "comp_library/IO_library/JSON/JSON_PARSER/what_ever.h"

namespace IO {
namespace JSONPARSER {
bool& operator<=(bool& lhs, const Whatever& w) {
  switch (w.type) {
    case WHATEVER_INTEGER:

      lhs = static_cast<bool>(w.whateverInteger);
      return lhs;

    case WHATEVER_BOOL:

      lhs = w.whateverBool;
      return lhs;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

int& operator<=(int& lhs, const Whatever& w) {
  switch (w.type) {
    case WHATEVER_INTEGER: {
      lhs = w.whateverInteger;
      return lhs;
    }

    default: { throw std::logic_error(__FUNCTION__); }
  }
}

float& operator<=(float& lhs, const Whatever& w) {
  switch (w.type) {
    case WHATEVER_INTEGER: {
      lhs = static_cast<float>(w.whateverInteger);
      return lhs;
    }

    case WHATEVER_DOUBLE: {
      lhs = static_cast<float>(w.whateverDouble);
      return lhs;
    }

    default: { throw std::logic_error(__FUNCTION__); }
  }
}

double& operator<=(double& lhs, const Whatever& w) {
  switch (w.type) {
    case WHATEVER_INTEGER: {
      lhs = static_cast<double>(w.whateverInteger);
      return lhs;
    }
    case WHATEVER_DOUBLE: {
      lhs = w.whateverDouble;
      return lhs;
    }

    default: { throw std::logic_error(__FUNCTION__); }
  }
}

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

std::string& operator<=(std::string& lhs, const Whatever& w) {
  switch (w.type) {
    case WHATEVER_STRING: {
      lhs = std::string(w.valueString.begin(), w.valueString.end());
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

}  // JSONPARSER
}  // IO

#endif  // COMP_LIBRARY_IO_LIBRARY_JSON_JSON_PARSER_JSON_OPERATORS_H
