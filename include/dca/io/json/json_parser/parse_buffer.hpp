// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// JSON parse buffer.

#ifndef DCA_IO_JSON_JSON_PARSER_PARSE_BUFFER_HPP
#define DCA_IO_JSON_JSON_PARSER_PARSE_BUFFER_HPP

#include <iostream>
#include <string>
#include <vector>

namespace dca {
namespace io {
// dca::io::

class ParseBuffer {
public:
  ParseBuffer() : theCharacters(), trace(false) {}

  void clear() {
    theCharacters.clear();

    if (trace)
      std::wcout << L"   ParseBuffer: clear()\n";
  }

  void put(wchar_t wc) {
    theCharacters.push_back(wc);
  }

  std::wstring str() {
    return std::wstring(theCharacters.begin(), theCharacters.end());
  }

  std::string to_string() {
    return std::string(theCharacters.begin(), theCharacters.end());
  }

public:
  std::vector<wchar_t> theCharacters;
  bool trace;
};

}  // io
}  // dca

#endif  // DCA_IO_JSON_JSON_PARSER_PARSE_BUFFER_HPP
