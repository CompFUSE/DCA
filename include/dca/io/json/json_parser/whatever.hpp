// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Whatever class.

#ifndef DCA_IO_JSON_JSON_PARSER_WHATEVER_HPP
#define DCA_IO_JSON_JSON_PARSER_WHATEVER_HPP

#include <cstdlib>
#include <cassert>
#include <map>
#include <sstream>
#include <vector>

#include "dca/io/json/json_parser/json_enumerations.hpp"

namespace dca {
namespace io {
// dca::io::

class Whatever {
public:
  typedef std::map<std::wstring, Whatever> WhateverMap;
  typedef std::vector<Whatever> WhateverVector;

  Whatever();

  Whatever(JSON_whatever_type t);

  std::size_t size() const {
    assert(type == WHATEVER_VECTOR);
    return whateverVector.size();
  }

  // std::wstring name()  const ;
  std::string name() const;

  Whatever& operator[](const std::wstring key) {
    assert(type == WHATEVER_MAP and whateverMap.count(key) == 1);
    return whateverMap[key];
  }
  Whatever& operator[](const std::string key) {
    std::wstring wKey(key.begin(), key.end());
    assert(type == WHATEVER_MAP and whateverMap.count(wKey) == 1);
    return whateverMap[wKey];
  }
  Whatever& operator[](std::size_t index) {
    assert(type == WHATEVER_VECTOR);
    assert(index < whateverVector.size());
    return whateverVector[index];
  }

  const Whatever& operator[](const std::wstring key) const {
    assert(type == WHATEVER_MAP and whateverMap.count(key) == 1);
    WhateverMap& wm = const_cast<WhateverMap&>(whateverMap);
    return wm[key];
  }
  const Whatever& operator[](const std::string key) const {
    std::wstring wKey(key.begin(), key.end());
    WhateverMap& wm = const_cast<WhateverMap&>(whateverMap);
    return wm[wKey];
  }
  const Whatever& operator[](std::size_t index) const {
    assert(type == WHATEVER_VECTOR);
    assert(index < whateverVector.size());
    return whateverVector[index];
  }

  static Whatever null();

  Whatever& back() {
    assert(type == WHATEVER_VECTOR);
    return whateverVector.back();
  }

  template <typename T>
  Whatever& push_back(T& value);
  Whatever& push_back();
  Whatever& push_back_null();

  template <typename T>
  Whatever& operator=(const T& value) {
    set(value);
    return *this;
  }

private:
  static std::wstring typeName(JSON_whatever_type t);
  static std::string ntypeName(JSON_whatever_type t);

  void collectName(std::wostringstream& nam) const;

  int count(std::wstring key) {
    assert(type == WHATEVER_VECTOR);
    return whateverMap.count(key);
  }

  template <typename T>
  void get(T& value) {
    value = valueString;
  }

  void setNull() {
    type = WHATEVER_NULL;
    valueString = L"NULL";
  }
  void set(const std::wstring& value) {
    type = WHATEVER_STRING;
    valueString = value;
  }
  void set(bool value) {
    type = WHATEVER_BOOL;
    whateverBool = value;
  }
  void set(int value) {
    type = WHATEVER_INTEGER;
    whateverInteger = value;
  }
  void set(double value) {
    type = WHATEVER_DOUBLE;
    whateverDouble = value;
  }

  // private:
public:
  JSON_whatever_type type;
  Whatever* parent;

  std::wstring valueString;

  WhateverMap whateverMap;
  WhateverVector whateverVector;
  bool whateverBool;
  int whateverInteger;
  double whateverDouble;

  std::wstring myKey;
  int myIndex;

  std::string filename;
  int startPos;
  int endPos;
};

template <typename T>
Whatever& Whatever::push_back(T& value) {
  assert(type == WHATEVER_VECTOR);

  whateverVector.push_back(Whatever());
  whateverVector.back() = value;

  return *this;
}

}  // io
}  // dca

#endif  // DCA_IO_JSON_JSON_PARSER_WHATEVER_HPP
