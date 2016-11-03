// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// JSON mode stack.

#ifndef DCA_IO_JSON_JSON_PARSER_JSON_MODE_STACK_HPP
#define DCA_IO_JSON_JSON_PARSER_JSON_MODE_STACK_HPP

#include <string>
#include <vector>

#include "dca/io/json/json_parser/json_enumerations.hpp"

namespace dca {
namespace io {
// dca::io::

class JSON_mode_stack {
public:
  JSON_mode_stack() : stack(1, MODE_DONE) {}

  std::string modeName(JSON_mode_type m);

  void push(JSON_mode_type mode) {
    stack.push_back(mode);
  }

  void pop(JSON_mode_type expectedMode);

  const JSON_mode_type& currentMode() {
    return stack.back();
  }

public:
  std::vector<JSON_mode_type> stack;
};

}  // io
}  // dca

#endif  //  DCA_IO_JSON_JSON_PARSER_JSON_MODE_STACK_HPP
