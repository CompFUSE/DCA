// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements json_operators.hpp.

#include "dca/io/json/json_parser/json_operators.hpp"

namespace dca {
namespace io {
// dca::io::

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

std::string& operator<=(std::string& lhs, const Whatever& w) {
  switch (w.type) {
    case WHATEVER_STRING: {
      lhs = std::string(w.valueString.begin(), w.valueString.end());
      return lhs;
    }

    default: { throw std::logic_error(__FUNCTION__); }
  }
}

}  // io
}  // dca
