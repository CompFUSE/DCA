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

#ifndef COMP_LIBRARY_IO_LIBRARY_JSON_JSON_PARSER_JSON_MODE_STACK_H
#define COMP_LIBRARY_IO_LIBRARY_JSON_JSON_PARSER_JSON_MODE_STACK_H

#include <sstream>
#include <string>
#include <vector>

#include "comp_library/IO_library/JSON/JSON_PARSER/JSON_enumerations.h"

namespace IO {
namespace JSONPARSER {
class JSON_mode_stack {
public:
  JSON_mode_stack();

  std::string modeName(JSON_mode_type m);

  void push(JSON_mode_type mode);

  void pop(JSON_mode_type expectedMode);

  const JSON_mode_type& currentMode();

public:
  std::vector<JSON_mode_type> stack;
};

JSON_mode_stack::JSON_mode_stack() : stack(1, MODE_DONE) {}

std::string JSON_mode_stack::modeName(JSON_mode_type m) {
  switch (m) {
    case MODE_ARRAY:
      return "MODE_ARRAY";
    case MODE_DONE:
      return "MODE_DONE";
    case MODE_KEY:
      return "MODE_KEY";
    case MODE_OBJECT:
      return "MODE_OBJECT";
    default:
      return "MODE_UNKNOWN";
  }
}

void JSON_mode_stack::push(JSON_mode_type mode) {
  stack.push_back(mode);
}

void JSON_mode_stack::pop(JSON_mode_type expectedMode) {
  if (stack.size() == 0) {
    std::ostringstream msg;
    msg << "JsonParser.pop was expecting mode " << modeName(expectedMode)
        << " to be on the back of the stack. \n"
        << "However the stack was empty!\n";
    throw std::logic_error(msg.str());
  }

  if (expectedMode != stack.back()) {
    std::ostringstream msg;
    msg << "JsonParser.pop was expecting mode " << modeName(expectedMode)
        << " to be on the back of the stack. \n"
        << "However the back of the stack contained " << modeName(stack.back()) << "\n";
    throw std::logic_error(msg.str());
  }

  stack.pop_back();
}

const JSON_mode_type& JSON_mode_stack::currentMode() {
  return stack.back();
}
}
}

#endif  // COMP_LIBRARY_IO_LIBRARY_JSON_JSON_PARSER_JSON_MODE_STACK_H
