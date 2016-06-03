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

#ifndef COMP_LIBRARY_IO_LIBRARY_JSON_JSON_PARSER_JSON_STATE_AND_ACTION_PAIR_H
#define COMP_LIBRARY_IO_LIBRARY_JSON_JSON_PARSER_JSON_STATE_AND_ACTION_PAIR_H

#include <vector>
#include "comp_library/IO_library/JSON/JSON_PARSER/JSON_enumerations.h"

namespace IO {
namespace JSONPARSER {
class state_and_action_pair {
public:
  state_and_action_pair(const state_and_action_pair& other);

  state_and_action_pair(JSON_state_type s, JSON_action_type a);

  state_and_action_pair(JSON_state_type s, JSON_action_type a, JSON_action_type b);

public:
  JSON_state_type newState;
  std::vector<JSON_action_type> actions;
};

state_and_action_pair::state_and_action_pair(const state_and_action_pair& other)
    : newState(other.newState), actions(other.actions) {}

state_and_action_pair::state_and_action_pair(JSON_state_type s, JSON_action_type a)
    : newState(s), actions(1) {
  actions[0] = a;
}

state_and_action_pair::state_and_action_pair(JSON_state_type s, JSON_action_type a, JSON_action_type b)
    : newState(s), actions(2) {
  actions[0] = a;
  actions[1] = b;
}
}
}

#endif  // COMP_LIBRARY_IO_LIBRARY_JSON_JSON_PARSER_JSON_STATE_AND_ACTION_PAIR_H
