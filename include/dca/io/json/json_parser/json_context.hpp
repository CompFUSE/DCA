// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// JSON context.

#ifndef DCA_IO_JSON_JSON_PARSER_JSON_CONTEXT_HPP
#define DCA_IO_JSON_JSON_PARSER_JSON_CONTEXT_HPP

#include <iostream>
#include <string>
#include <vector>

#include "dca/io/json/json_parser/json_operators.hpp"
#include "dca/io/json/json_parser/parse_buffer.hpp"
#include "dca/io/json/json_parser/whatever.hpp"

namespace dca {
namespace io {
// dca::io::

class JSON_context {
public:
  JSON_context() : result(), stack(1, &result), key(L""), trace(false) {}

  void begin_object() {
    beginObjectOrArray<WHATEVER_MAP>();
  }
  void end_object() {
    endObjectOrArray<WHATEVER_MAP>();
  }

  void begin_numeric_array(std::string filename, size_t charNum) {
    beginObjectOrArray<WHATEVER_MAT>();
    currentObject().filename = filename;
    currentObject().startPos = charNum;
  }
  void end_numeric_array(size_t charNum) {
    currentObject().endPos = charNum;
    endObjectOrArray<WHATEVER_MAT>();
  }

  void begin_array() {
    beginObjectOrArray<WHATEVER_VECTOR>();
  }
  void end_array() {
    endObjectOrArray<WHATEVER_VECTOR>();
  }

  void Integer(ParseBuffer& s) {
    int i;
    s >> i;
    set(i);
  }
  void Float(ParseBuffer& s) {
    double d;
    s >> d;
    set(d);
  }
  void Null() {
    set(Whatever::null());
  }
  void True() {
    set(true);
  }
  void False() {
    set(false);
  }
  void String(ParseBuffer& s) {
    set(s.str());
  }
  void Key(const std::wstring& s) {
    key = s;
    if (trace)
      std::wcout << L"   key = '" << key << L"'\n";
  }

private:
  Whatever& currentObject() {
    return *stack.back();
  }

  // Object and Array Referenced objects are created when needed on LHS of assignment!
  Whatever& referencedObject();

  template <JSON_whatever_type type_id>
  void beginObjectOrArray();

  template <int TypeId>
  void endObjectOrArray();

  std::string CurrentContext() {
    return referencedObject().name();
  }

  template <typename T>
  void set(const T& v) {
    Whatever& refObj(referencedObject());
    refObj = v;
  }

public:
  Whatever result;
  std::vector<Whatever*> stack;
  std::wstring key;
  bool trace;
};

template <JSON_whatever_type TypeId>
void JSON_context::beginObjectOrArray() {
  Whatever& refObj = referencedObject();  // Generally creates the object or array
                                          // (except for the first time when refObject == result )
  refObj.type = TypeId;                   // Sets the type
  if (&refObj != &result)      // result is already on the stack, so we dont' need to push it on.
    stack.push_back(&refObj);  // Make the referenced object the current object

  if (trace)
    std::cout << " Set the type of " << refObj.name() << " to " << name(TypeId)
              << " and make it the current object\n";
}

template <int TypeId>
void JSON_context::endObjectOrArray() {
  if (trace)
    std::cout << "   JSON_context is ending object '" << currentObject().name()
              << "' by popping the object stack!\n";

  assert(currentObject().type == TypeId);

  stack.pop_back();

  if (trace) {
    if (stack.size() > 0) {
      std::cout << "   current object is now '" << currentObject().name() << "\n";
    }
    else
      std::cout << "   Parsing completed! \n ";
  }
}

}  // io
}  // dca

#endif  // DCA_IO_JSON_JSON_PARSER_JSON_CONTEXT_HPP
