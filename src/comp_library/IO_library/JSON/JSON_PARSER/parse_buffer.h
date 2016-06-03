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

#ifndef COMP_LIBRARY_IO_LIBRARY_JSON_JSON_PARSER_PARSE_BUFFER_H
#define COMP_LIBRARY_IO_LIBRARY_JSON_JSON_PARSER_PARSE_BUFFER_H

#include <iostream>
#include <string>
#include <vector>

namespace IO {
namespace JSONPARSER {
class ParseBuffer {
public:
  ParseBuffer();

  void clear();

  void put(wchar_t wc);

  std::wstring str();

  std::string to_string();

public:
  std::vector<wchar_t> theCharacters;
  bool trace;
};

ParseBuffer::ParseBuffer() : theCharacters(), trace(false) {}

void ParseBuffer::clear() {
  theCharacters.clear();

  if (trace)
    std::wcout << L"   ParseBuffer: clear()\n";
}

void ParseBuffer::put(wchar_t wc) {
  theCharacters.push_back(wc);
}

std::wstring ParseBuffer::str() {
  return std::wstring(theCharacters.begin(), theCharacters.end());
}

std::string ParseBuffer::to_string() {
  return std::string(theCharacters.begin(), theCharacters.end());
}
}
}

#endif  // COMP_LIBRARY_IO_LIBRARY_JSON_JSON_PARSER_PARSE_BUFFER_H
