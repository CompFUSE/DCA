// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef DCA_IO_JSON_JSON_PARSER_JSON_TRANSLATION_TABLE_HPP
#define DCA_IO_JSON_JSON_PARSER_JSON_TRANSLATION_TABLE_HPP

#include <cassert>
#include <stdexcept>

#include "dca/io/json/json_parser/json_character_mapper.hpp"
#include "dca/io/json/json_parser/json_enumerations.hpp"
#include "dca/io/json/json_parser/state_and_action_pair.hpp"

namespace dca {
namespace io {
// dca::io::

class JSON_translation_table : public JSON_character_mapper {
public:
  static state_and_action_pair get_state_and_action_pair(JSON_state_type& state,
                                                         JSON_character_class_type& cls);

private:
  static bool check_state_and_class(JSON_state_type& state, JSON_character_class_type& cls);

  static state_and_action_pair begin_object_or_array(JSON_character_class_type& cls);

  static state_and_action_pair begin_value(JSON_character_class_type& cls);
  static state_and_action_pair end_value(JSON_character_class_type& cls);

  static state_and_action_pair begin_string(JSON_character_class_type& cls);
  static state_and_action_pair read_string(JSON_character_class_type& cls);

  static state_and_action_pair read_integer(JSON_character_class_type& cls);
  static state_and_action_pair read_fraction(JSON_character_class_type& cls);

  static state_and_action_pair begin_exponent(JSON_character_class_type& cls);
  static state_and_action_pair read_exponent(JSON_character_class_type& cls);

  static state_and_action_pair read_NULL(JSON_state_type& state, JSON_character_class_type& cls);

  static state_and_action_pair read_false(JSON_state_type& state, JSON_character_class_type& cls);
  static state_and_action_pair read_true(JSON_state_type& state, JSON_character_class_type& cls);

  static state_and_action_pair read_array(JSON_state_type& state, JSON_character_class_type& cls);
};

}  // io
}  // dca

#endif  // DCA_IO_JSON_JSON_PARSER_JSON_TRANSLATION_TABLE_HPP
