// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements json_enumerations.hpp.

#include "dca/io/json/json_parser/json_enumerations.hpp"

namespace dca {
namespace io {
// dca::io::

std::string name(JSON_character_class_type cls) {
  switch (cls) {
    case C_SPACE:
      return "C_SPACE space ";
    case C_WHITE:
      return "C_WHITE other whitespace ";
    case C_LPARN:
      return "C_LPARN ( ";
    case C_RPARN:
      return "C_RPARN ) ";
    case C_LCURB:
      return "C_LCURB { ";
    case C_RCURB:
      return "C_RCURB } ";
    case C_LSQRB:
      return "C_LSQRB [ ";
    case C_RSQRB:
      return "C_RSQRB ] ";
    case C_COLON:
      return "C_COLON : ";
    case C_COMMA:
      return "C_COMMA , ";
    case C_QUOTE:
      return "C_QUOTE \\\" ";
    case C_BACKS:
      return "C_BACKS \\ ";
    case C_SLASH:
      return "C_SLASH / ";
    case C_PLUS:
      return "C_PLUS + ";
    case C_MINUS:
      return "C_MINUS - ";
    case C_POINT:
      return "C_POINT . ";
    case C_ZERO:
      return "C_ZERO 0 ";
    case C_DIGIT:
      return "C_DIGIT 123456789 ";
    case C_LOW_A:
      return "C_LOW_A a ";
    case C_LOW_B:
      return "C_LOW_B b ";
    case C_LOW_C:
      return "C_LOW_C c ";
    case C_LOW_D:
      return "C_LOW_D d ";
    case C_LOW_E:
      return "C_LOW_E e ";
    case C_LOW_F:
      return "C_LOW_F f ";
    case C_LOW_L:
      return "C_LOW_L l ";
    case C_LOW_N:
      return "C_LOW_N n ";
    case C_LOW_R:
      return "C_LOW_R r ";
    case C_LOW_S:
      return "C_LOW_S s ";
    case C_LOW_T:
      return "C_LOW_T t ";
    case C_LOW_U:
      return "C_LOW_U u ";
    case C_LOW_Y:
      return "C_LOW_Y y ";
    case C_ABCDF:
      return "C_ABCDF ABCDF ";
    case C_E:
      return "C_E E ";
    case C_ETC:
      return "C_ETC everything else ";
    case C_STAR:
      return "C_STAR * ";
    case C_EOF:
      return "C_ERR error ";
    case C_ERR:
      return "C_ERR error ";
    case NR_CLASSES:
      return "NP_CLASSES ";
    default:
      return "Unkown CharacterClass ";
  }
}

std::string name(JSON_action_type& action) {
  switch (action) {
    case Consume:
      return " Consume the character    ";
    case BeginObject:
      return " BeginObject              ";
    case BeginArray:
      return " BeginArray               ";
    case BeginMatrix:
      return " BeginMatrix               ";
    case EndObject:
      return " endObject              ";
    case EndArray:
      return " EndArray               ";
    case EndMatrix:
      return " EndMatrix               ";
    case DoNext:
      return " DoNext                 ";
    case RecordKey:
      return " RecordKey                ";
    case RecordChar:
      return " RecordChar               ";
    case RecordString:
      return " RecordString             ";
    case RecordInteger:
      return " RecordInteger            ";
    case RecordFloat:
      return " RecordFloat            ";
    case RecordTrue:
      return " RecordTrue            ";
    case RecordFalse:
      return " RecordFalse            ";
    case RecordNull:
      return " RecordNull            ";
    case EndFile:
      return " EndFile            ";
    case Abort:
      return " Abort Parsing            ";
    default:
      return " Unkown action code";
  }
}

std::string name(JSON_state_type& state) {
  switch (state) {
    case GO:
      return "GO: start    ";
    case VA:
      return "VA: looking for a value    ";
    case EV:
      return "EV: Looking after a value  ";
    case BS:
      return "BS: Looking for the begining of a string  ";
    case ST:
      return "ST: Looking for string characters         ";
    case IT:
      return "IT: Looking for an integer  ";
    case FR:
      return "FR: Looking for the integer after the point ";
    case EX:
      return "EX: Looking for the exponent ";
    case EX2:
      return "EX2: Looking for the integer part of the exponent ";
    case T:
      return "T: Looking for _rue (true)   ";
    case TR:
      return "TR: Looking for __ue (true)        ";
    case TRU:
      return "TRU: Looking for ___e (true)       ";
    case F:
      return "F: Looking for _alse (false)   ";
    case FA:
      return "FA: Looking for __lse (false)   ";
    case FAL:
      return "FAL: Looking for ___se (false)   ";
    case FALS:
      return "FALS: Looking for ____e (false)   ";
    case N:
      return "N: Looking for _ull (null)     ";
    case NU:
      return "NU: Looking for __ll (null)     ";
    case NUL:
      return "NUL: Looking for ___l (null)     ";
    case A:
      return "A: Looking for _rray (array)     ";
    case AR:
      return "AR: Looking for __ray (array)     ";
    case ARR:
      return "ARR: Looking for ___ay (array)     ";
    case ARRA:
      return "ARRA: Looking for ___y (array)     ";
    case ARRAY:
      return "ARRAY: Looking for ____( (array)     ";
    case MS:
      return "skipping matrix chars     ";
    case END:
      return "END: End of Parsing                ";
    default:
      return " Unkown state code";
  }
}

std::string name(JSON_whatever_type t) {
  switch (t) {
    case WHATEVER_MAT:
      return "WHATEVER_MAT";
    case WHATEVER_MAP:
      return "WHATEVER_MAP";
    case WHATEVER_VECTOR:
      return "WHATEVER_VECTOR";
    case WHATEVER_MATRIX:
      return "WHATEVER_MATRIX";
    case WHATEVER_STRING:
      return "WHATEVER_STRING";
    case WHATEVER_INTEGER:
      return "WHATEVER_INTEGER";
    case WHATEVER_DOUBLE:
      return "WHATEVER_DOUBLE";
    case WHATEVER_BOOL:
      return "WHATEVER_BOOL";
    case WHATEVER_UNKNOWN:
      return "WHATEVER_UNKNOWN";
    default:
      throw std::logic_error("Whatever given wrong type");
  }
}

}  // io
}  // dca
