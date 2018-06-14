// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements whatever.hpp.

#include "dca/io/json/json_parser/whatever.hpp"

namespace dca {
namespace io {
// dca::io::

Whatever::Whatever()
    : type(WHATEVER_UNKNOWN),
      parent(0),
      valueString(),
      whateverMap(),
      whateverVector(),
      whateverBool(true),
      whateverInteger(0),
      whateverDouble(0),
      myKey(L"?"),
      myIndex(-1),
      startPos(0),
      endPos(0) {}

Whatever::Whatever(JSON_whatever_type t)
    : type(t),
      parent(0),
      valueString(),
      whateverMap(),
      whateverVector(),
      whateverBool(true),
      whateverInteger(0),
      whateverDouble(0),
      myKey(L"?"),
      myIndex(-1),
      startPos(0),
      endPos(0) {}

std::string Whatever::name() const {
  std::wostringstream nam;

  collectName(nam);
  nam << "{" << typeName(type) << "}";

  std::wstring wname = nam.str();  //(name());

  return std::string(wname.begin(), wname.end());
}

Whatever Whatever::null() {
  Whatever result;
  result.setNull();
  return result;
}

Whatever& Whatever::push_back() {
  assert(type == WHATEVER_VECTOR);
  whateverVector.push_back(Whatever());
  return *this;
}

Whatever& Whatever::push_back_null() {
  assert(type == WHATEVER_VECTOR);

  whateverVector.push_back(Whatever());
  whateverVector.back().setNull();

  return *this;
}

std::wstring Whatever::typeName(JSON_whatever_type t) {
  switch (t) {
    case WHATEVER_MAT:
      return L"WHATEVER_MAT";
    case WHATEVER_MAP:
      return L"WHATEVER_MAP";
    case WHATEVER_VECTOR:
      return L"WHATEVER_VECTOR";
    case WHATEVER_MATRIX:
      return L"WHATEVER_MATRIX";
    case WHATEVER_STRING:
      return L"WHATEVER_STRING";
    case WHATEVER_INTEGER:
      return L"WHATEVER_INTEGER";
    case WHATEVER_DOUBLE:
      return L"WHATEVER_DOUBLE";
    case WHATEVER_BOOL:
      return L"WHATEVER_BOOL";
    case WHATEVER_UNKNOWN:
      return L"WHATEVER_UNKNOWN";
    default:
      throw std::logic_error("Whatever::typeName given wrong type");
  }
}

std::string Whatever::ntypeName(JSON_whatever_type t) {
  std::wstring wname(typeName(t));
  return std::string(wname.begin(), wname.end());
}

void Whatever::collectName(std::wostringstream& nam) const {
  if (parent == 0) {
    nam << "root";
    return;
  }

  parent->collectName(nam);

  switch (parent->type) {
    case WHATEVER_MAP:
      nam << L"[" << myKey << L"]";
      break;
    case WHATEVER_VECTOR:
      nam << L"[" << myIndex << L"]";
      break;
    default:
      nam << L"[ ?? not vector or map! ]";
  }
}

}  // io
}  // dca
