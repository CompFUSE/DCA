// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// JSON parser.

#ifndef DCA_IO_JSON_JSON_PARSER_JSON_PARSER_HPP
#define DCA_IO_JSON_JSON_PARSER_JSON_PARSER_HPP

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "dca/io/json/json_parser/json_character_mapper.hpp"
#include "dca/io/json/json_parser/json_mode_stack.hpp"
#include "dca/io/json/json_parser/json_translation_table.hpp"
#include "dca/io/json/json_parser/parse_buffer.hpp"
#include "dca/io/json/json_parser/state_and_action_pair.hpp"

namespace dca {
namespace io {
// dca::io::

template <typename context_type>
class JSON_parser : public JSON_character_mapper, public JSON_mode_stack {
public:
  JSON_parser();

  bool execute(std::wistream& inputStream);

  context_type& get_JSON_tree();

private:
  std::pair<wchar_t, JSON_character_class_type> get_next_character_and_class(std::wistream& inputStream);

  void begin_entity(JSON_action_type& action);
  void end_entity(JSON_action_type& action);

  void record_entity(JSON_action_type& action);

  void abort_parsing(wchar_t nextChar, JSON_character_class_type nextClass, JSON_state_type state,
                     JSON_action_type action);

  void setCommentState();

  void unloadBuffer();

  void performAction(wchar_t nextChar, JSON_character_class_type nextClass, JSON_state_type state,
                     JSON_action_type action);

private:
  context_type ctx;
  ParseBuffer buffer;
  JSON_state_type state;

  JSON_value_type current_value_type;

  bool allow_comments;
  JSON_state_type before_comment_state;
  bool comment;
  size_t comment_begin_offset;

  size_t numChar;
  size_t numLines;

  bool trace;

  std::string filename;
};

template <typename context_type>
JSON_parser<context_type>::JSON_parser()
    : ctx(),
      buffer(),
      state(GO),
      current_value_type(JSON_T_NONE),
      // escaped              (false),
      allow_comments(true),
      before_comment_state(GO),
      comment(false),
      comment_begin_offset(0),
      numChar(0),
      numLines(0),
      trace(false),
      filename("") {}

template <typename context_type>
context_type& JSON_parser<context_type>::get_JSON_tree() {
  return ctx;
}

template <typename context_type>
void JSON_parser<context_type>::setCommentState() {
  switch (currentMode()) {
    case MODE_ARRAY:
    case MODE_OBJECT:
      switch (state) {
        case VA:
        case AR:
          before_comment_state = state;
          break;
        default:
          before_comment_state = GO;
          break;
      }
      break;

    default:
      before_comment_state = state;
      break;
  }

  comment = true;
}

template <typename context_type>
void JSON_parser<context_type>::unloadBuffer() {
  switch (current_value_type) {
    case JSON_T_STRING:
      ctx.String(buffer);
      break;
    default:
      break;
  }

  current_value_type = JSON_T_NONE;

  buffer.clear();
}

template <typename context_type>
void JSON_parser<context_type>::performAction(wchar_t nextChar, JSON_character_class_type nextClass,
                                              JSON_state_type state, JSON_action_type action) {
  switch (action) {
    case Consume:
      return;

    case DoNext: {
      unloadBuffer();
      return;
    }

    case RecordChar: {
      buffer.put(nextChar);
      return;
    }

    case BeginArray:
    case BeginObject:
    case BeginMatrix: {
      begin_entity(action);
      return;
    }

    case EndArray:
    case EndObject:
    case EndMatrix: {
      end_entity(action);
      return;
    }

    case RecordString:
    case RecordKey:
    case RecordFloat:
    case RecordInteger:
    case RecordTrue:
    case RecordFalse:
    case RecordNull: {
      record_entity(action);
      break;
    }

    case Abort: {
      abort_parsing(nextChar, nextClass, state, action);
      break;
    }

    case EndFile: {
      std::cout << "\n\n\t Parsing completed! read " << numChar << " characters and " << numLines
                << " lines.\n";
      break;
    }

    default: { throw std::logic_error(__FUNCTION__); }
  }
}

template <typename context_type>
void JSON_parser<context_type>::begin_entity(JSON_action_type& action) {
  switch (action) {
    case BeginObject:
      ctx.begin_object();
      break;

    case BeginArray:
      ctx.begin_array();
      break;

    case BeginMatrix:
      ctx.begin_numeric_array(filename, numChar);
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }

  buffer.clear();
}

template <typename context_type>
void JSON_parser<context_type>::end_entity(JSON_action_type& action) {
  switch (action) {
    case EndObject:
      unloadBuffer();
      ctx.end_object();
      break;

    case EndArray:
      unloadBuffer();
      ctx.end_array();
      break;

    case EndMatrix:
      ctx.end_numeric_array(numChar);
      buffer.clear();
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

template <typename context_type>
void JSON_parser<context_type>::record_entity(JSON_action_type& action) {
  switch (action) {
    case RecordString:
      current_value_type = JSON_T_STRING;
      break;

    case RecordKey:
      current_value_type = JSON_T_NONE;
      ctx.Key(buffer.str());
      buffer.clear();
      break;

    case RecordFloat:
      current_value_type = JSON_T_FLOAT;
      ctx.Float(buffer);
      buffer.clear();
      break;

    case RecordInteger:
      current_value_type = JSON_T_INTEGER;
      ctx.Integer(buffer);
      buffer.clear();
      break;

    case RecordTrue:
      current_value_type = JSON_T_TRUE;
      ctx.True();
      break;

    case RecordFalse:
      current_value_type = JSON_T_FALSE;
      ctx.False();
      break;

    case RecordNull:
      current_value_type = JSON_T_NULL;
      ctx.Null();
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

template <typename context_type>
void JSON_parser<context_type>::abort_parsing(wchar_t nextChar, JSON_character_class_type nextClass,
                                              JSON_state_type state, JSON_action_type action) {
  std::cout << "JsonParser::performAction was sent abort from JSON_parser : \n"
            << "  nextChar   = '" << ((char)wctob(nextChar)) << "'\n"
            << "  nextClass  = " << name(nextClass) << "\n"
            << "  state      = " << name(state) << "\n"
            << "  action     = " << name(action) << "\n"
            << "  character# = " << numChar << "\n"
            << "  line#      = " << numLines << "\n";

  std::cout << "\n\t the error was recorded in the buffer : " << buffer.to_string() << "\n\n";
  throw std::logic_error(__FUNCTION__);
}

template <typename context_type>
bool JSON_parser<context_type>::execute(std::wistream& inputStream) {
  std::pair<wchar_t, JSON_character_class_type> next = get_next_character_and_class(inputStream);

  wchar_t& nextChar(next.first);
  JSON_character_class_type& nextClass(next.second);

  state_and_action_pair pair = JSON_translation_table::get_state_and_action_pair(state, nextClass);

  std::vector<JSON_action_type>& actions = pair.actions;

  if (trace)
    std::cout << "actions.size() = " << actions.size();

  for (size_t i = 0; i < actions.size(); i++) {
    performAction(nextChar, nextClass, state, actions[i]);

    if (actions[i] == EndFile)
      return false;
  }

  state = pair.newState;

  return !inputStream.eof();
}

template <typename context_type>
std::pair<wchar_t, JSON_character_class_type> JSON_parser<context_type>::get_next_character_and_class(
    std::wistream& inputStream) {
  std::pair<wchar_t, JSON_character_class_type> result;

  wchar_t& nextChar = result.first;
  JSON_character_class_type& nextClass = result.second;

  do {
    if (inputStream.eof()) {
      nextChar = EOF;
      nextClass = C_EOF;

      return result;
    }

    nextChar = inputStream.get();
    if (inputStream.eof()) {
      nextChar = EOF;
      nextClass = C_EOF;

      return result;
    }
    numChar++;

    if (nextChar == L'\n')
      numLines++;

    if (nextChar == L'\\' and state == ST) {
      nextChar = get_escaped_character(inputStream);
      numChar++;
    }

    nextClass = map_char_to_class(nextChar);

    if (state == ST || comment)
      return result;

    if (nextChar < 0)
      nextClass = C_WHITE;
  } while (is_white_space(nextClass));

  assert(nextChar >= 0);

  return result;
}

}  // io
}  // dca

#endif  // DCA_IO_JSON_JSON_PARSER_JSON_PARSER_HPP
