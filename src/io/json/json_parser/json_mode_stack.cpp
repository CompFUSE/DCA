// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements json_mode_stack.hpp.

#include "dca/io/json/json_parser/json_mode_stack.hpp"
#include <sstream>

namespace dca {
namespace io {
// dca::io::

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

}  // io
}  // dca
