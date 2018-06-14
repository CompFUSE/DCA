// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements json_character_mapper.hpp.

#include "dca/io/json/json_parser/json_character_mapper.hpp"

namespace dca {
namespace io {
// dca::io::

// This array maps the 128 ASCII characters into character classes.
// The remaining Unicode characters should be mapped to C_ETC.
// Non-whitespace control characters are errors.
JSON_character_class_type JSON_character_mapper::ascii_class[] = {
    C_ERR,   C_ERR,   C_ERR,   C_ERR,   C_ERR,   C_ERR,   C_ERR,   C_ERR,
    C_ERR,   C_WHITE, C_WHITE, C_ERR,   C_ERR,   C_WHITE, C_ERR,   C_ERR,
    C_ERR,   C_ERR,   C_ERR,   C_ERR,   C_ERR,   C_ERR,   C_ERR,   C_ERR,
    C_ERR,   C_ERR,   C_ERR,   C_ERR,   C_ERR,   C_ERR,   C_ERR,   C_ERR,

    C_SPACE, C_ETC,   C_QUOTE, C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,
    C_LPARN, C_RPARN, C_STAR,  C_PLUS,  C_COMMA, C_MINUS, C_POINT, C_SLASH,
    C_ZERO,  C_DIGIT, C_DIGIT, C_DIGIT, C_DIGIT, C_DIGIT, C_DIGIT, C_DIGIT,
    C_DIGIT, C_DIGIT, C_COLON, C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,

    C_ETC,   C_ABCDF, C_ABCDF, C_ABCDF, C_ABCDF, C_E,     C_ABCDF, C_ETC,
    C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,
    C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,
    C_ETC,   C_ETC,   C_ETC,   C_LSQRB, C_BACKS, C_RSQRB, C_ETC,   C_ETC,

    C_ETC,   C_LOW_A, C_LOW_B, C_LOW_C, C_LOW_D, C_LOW_E, C_LOW_F, C_ETC,
    C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_LOW_L, C_ETC,   C_LOW_N, C_ETC,
    C_ETC,   C_ETC,   C_LOW_R, C_LOW_S, C_LOW_T, C_LOW_U, C_ETC,   C_ETC,
    C_ETC,   C_LOW_Y, C_ETC,   C_LCURB, C_ETC,   C_RCURB, C_ETC,   C_ETC};

JSON_character_class_type JSON_character_mapper::map_char_to_class(wchar_t widec) {
  if (widec >= 128)
    return C_ETC;
  else {
    JSON_character_class_type result = ascii_class[wctob(widec)];

    switch (result) {
      case C_ERR: {
        std::cout << "\n\n\t" << L"JSON_character_mapper::wcharToClass(" << widec
                  << L") resulted in an error!\n\n\n";
        throw std::logic_error(__FUNCTION__);
      }

      default:

        return result;
    }
  }
}

bool JSON_character_mapper::is_white_space(JSON_character_class_type& nextClass) {
  switch (nextClass) {
    case C_SPACE:
    case C_WHITE:
      return true;

    default:
      return false;
  }
}

wchar_t JSON_character_mapper::get_escaped_character(std::wistream& inputStream) {
  wchar_t secondChar = inputStream.get();

  switch (secondChar) {
    case L'b':
      return L'\b';
    case L'f':
      return L'\f';
    case L'n':
      return L'\n';
    case L'r':
      return L'\r';
    case L't':
      return L'\t';
    case L'"':
      return L'"';
    case L'\\':
      return L'\\';
    case L'/':
      return L'/';
    default:
      throw std::logic_error(
          "JsonParser: Encountered an escapped character that was not recognized.");
  }
}

}  // io
}  // dca
