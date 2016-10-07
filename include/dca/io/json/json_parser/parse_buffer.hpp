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
