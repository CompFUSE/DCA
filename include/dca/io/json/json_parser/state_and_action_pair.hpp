// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// State and action pair for JSON parser.

#ifndef DCA_IO_JSON_JSON_PARSER_STATE_AND_ACTION_PAIR_HPP
#define DCA_IO_JSON_JSON_PARSER_STATE_AND_ACTION_PAIR_HPP

#include <vector>
#include "dca/io/json/json_parser/json_enumerations.hpp"

namespace dca {
namespace io {
// dca::io::

class state_and_action_pair {
public:
  state_and_action_pair(const state_and_action_pair& other)
      : newState(other.newState), actions(other.actions) {}

  state_and_action_pair(JSON_state_type s, JSON_action_type a) : newState(s), actions(1) {
    actions[0] = a;
  }

  state_and_action_pair(JSON_state_type s, JSON_action_type a, JSON_action_type b)
      : newState(s), actions(2) {
    actions[0] = a;
    actions[1] = b;
  }

public:
  JSON_state_type newState;
  std::vector<JSON_action_type> actions;
};

}  // io
}  // dca

#endif  // DCA_IO_JSON_JSON_PARSER_STATE_AND_ACTION_PAIR_HPP
