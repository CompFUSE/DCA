// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides various enumerations for the JSON parser.

#ifndef DCA_IO_JSON_JSON_PARSER_JSON_ENUMERATIONS_HPP
#define DCA_IO_JSON_JSON_PARSER_JSON_ENUMERATIONS_HPP

#include <stdexcept>
#include <string>

namespace dca {
namespace io {
// dca::io::

//  Characters are mapped into these 31 character classes. This allows for a significant reduction
//  in the size of the state transition table.
enum JSON_character_class {
  C_SPACE, /* space */
  C_WHITE, /* other whitespace */
  C_LPARN, /* (  */
  C_RPARN, /* ) */
  C_LCURB, /* {  */
  C_RCURB, /* } */
  C_LSQRB, /* [ */
  C_RSQRB, /* ] */
  C_COLON, /* : */
  C_COMMA, /* , */
  C_QUOTE, /* " */
  C_BACKS, /* \ */
  C_SLASH, /* / */
  C_PLUS,  /* + */
  C_MINUS, /* - */
  C_POINT, /* . */
  C_ZERO,  /* 0 */
  C_DIGIT, /* 123456789 */
  C_LOW_A, /* a */
  C_LOW_B, /* b */
  C_LOW_C, /* c */
  C_LOW_D, /* d */
  C_LOW_E, /* e */
  C_LOW_F, /* f */
  C_LOW_L, /* l */
  C_LOW_N, /* n */
  C_LOW_R, /* r */
  C_LOW_S, /* s */
  C_LOW_T, /* t */
  C_LOW_U, /* u */
  C_LOW_Y, /* y */
  C_ABCDF, /* ABCDF */
  C_E,     /* E */
  C_ETC,   /* everything else */
  C_STAR,  /* * */
  C_EOF,   /* end of file */
  C_ERR,   /* error */
  NR_CLASSES
};

enum JSON_action {
  Consume,       /* Consume the character    */
  BeginObject,   /* BeginObject              */
  BeginArray,    /* BeginArray               */
  BeginMatrix,   /* BeginMatrix              */
  EndObject,     /* endObject              */
  EndArray,      /* EndArray               */
  EndMatrix,     /* EndMatrix               */
  DoNext,        /* DoNext                 */
  RecordKey,     /* RecordKey                */
  RecordChar,    /* RecordChar               */
  RecordString,  /* RecordString             */
  RecordInteger, /* RecordInteger            */
  RecordFloat,   /* RecordFloat            */
  RecordTrue,    /* RecordTrue            */
  RecordFalse,   /* RecordFalse            */
  RecordNull,    /* RecordNull            */
  EndFile,       /* EndFile               */
  Abort,         /* Abort Parsing            */
  NR_ACTIONS
};

enum JSON_mode { MODE_ARRAY = 1, MODE_DONE = 2, MODE_KEY = 3, MODE_OBJECT = 4 };

enum JSON_state {
  GO, /* start    */
  VA, /* looking for a value    */
  EV, /* Looking after a value  */
  BS, /* Looking for the begining of a string  */
  ST, /* Looking for string characters         */
  IT, /* Looking for an integer  */
  FR, /* Looking for the integer after the point */

  EX,  /* Looking for the exponent */
  EX2, /* Looking for the integer part of the exponent */

  T,   /* Looking for _rue (true)   */
  TR,  /* Looking for __ue (true)        */
  TRU, /* Looking for ___e (true)       */

  F,    /* Looking for _alse (false)   */
  FA,   /* Looking for __lse (false)   */
  FAL,  /* Looking for ___se (false)   */
  FALS, /* Looking for ____e (false)   */

  N,   /* Looking for _ull (null)     */
  NU,  /* Looking for __ll (null)     */
  NUL, /* Looking for ___l (null)     */

  A,     /* Looking for _rray (array)     */
  AR,    /* Looking for __ray (array)     */
  ARR,   /* Looking for ___ay (array)     */
  ARRA,  /* Looking for ___y (array)     */
  ARRAY, /* Looking for ____( (array)     */

  MS, /* Skipping matrix chars     */

  END, /* END of Parsing               */
  NR_STATES
};

enum JSON_value {
  JSON_T_NONE = 0,
  JSON_T_ARRAY_BEGIN,
  JSON_T_ARRAY_END,
  JSON_T_OBJECT_BEGIN,
  JSON_T_OBJECT_END,
  JSON_T_INTEGER,
  JSON_T_FLOAT,
  JSON_T_NULL,
  JSON_T_TRUE,
  JSON_T_FALSE,
  JSON_T_STRING,
  JSON_T_KEY,
  JSON_T_MAX
};

enum JSON_whatever {
  WHATEVER_MAT,
  WHATEVER_MAP,
  WHATEVER_VECTOR,
  WHATEVER_MATRIX,
  WHATEVER_STRING,
  WHATEVER_INTEGER,
  WHATEVER_DOUBLE,
  WHATEVER_BOOL,
  WHATEVER_NULL,
  WHATEVER_UNKNOWN
};

typedef JSON_mode JSON_mode_type;
typedef JSON_state JSON_state_type;
typedef JSON_value JSON_value_type;
typedef JSON_action JSON_action_type;
typedef JSON_character_class JSON_character_class_type;
typedef JSON_whatever JSON_whatever_type;

std::string name(JSON_character_class_type cls);
std::string name(JSON_action_type& action);
std::string name(JSON_state_type& state);
std::string name(JSON_whatever_type t);

}  // io
}  // dca

#endif  // DCA_IO_JSON_JSON_PARSER_JSON_ENUMERATIONS_HPP
